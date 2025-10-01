import instaseis
from tqdm import tqdm
import h5py
import os
from obspy import read, Stream, Trace
from obspy.imaging.beachball import mt2plane, MomentTensor
from cmt import compute_seismic_moment
import numpy as np
import json
import pyproj
from collections import defaultdict
from scipy.signal import hann


def fft_reconvolve_stf(db, trace_data, new_stf):
    """
    Deconvolve the original STF from a seismogram and convolve a new STF using FFT.

    Parameters:
    ----------
    db : instaseis.Database
        Instaseis database object, must have .info.sliprate and .info.dt.
    trace_data : np.ndarray
        Time series data (1D array) from Instaseis.
    new_stf : np.ndarray
        New STF to convolve with, e.g. moment-rate function. Must be bandlimited.
    Returns:
    --------
    np.ndarray
        Modified seismogram after deconvolution and convolution.
    """

    original_stf = db.info.sliprate
    dt = db.info.dt
    L = len(trace_data)
    nfft = int(2 ** np.ceil(np.log2(L + len(new_stf) - 1)))

    # FFT of original STF (with zero-padding)
    orig_stf_f = np.fft.rfft(original_stf, n=nfft)
    imax = np.argmax(original_stf)
    imin = np.argmin(original_stf)
    time_shift = 0.5 * (imax + imin) * dt
    if time_shift != 0.0:
        freqs = np.fft.rfftfreq(nfft, dt)
        orig_stf_f *= np.exp(1j * 2 * np.pi * freqs * time_shift)

    # FFT of new STF (with zero-padding)
    new_stf_f = np.fft.rfft(new_stf, n=nfft)

    # Apply taper to trace data to reduce ringing
    taper = np.ones_like(trace_data)
    taper_len = max(int(0.05 * L), 5)
    taper[-taper_len:] = hann(2 * taper_len)[taper_len:]
    tapered_trace = trace_data * taper

    # FFT of trace data
    data_f = np.fft.rfft(tapered_trace, n=nfft)

    # Deconvolve: avoid division by zero
    mask = np.abs(orig_stf_f) > 1e-12
    recon_f = np.zeros_like(data_f)
    recon_f[mask] = data_f[mask] / orig_stf_f[mask]
    recon_f[~mask] = 0.0

    # Convolve with new STF in frequency domain
    final_f = recon_f * new_stf_f

    # Inverse FFT and truncate to original length
    result = np.fft.irfft(final_f)[:L]
    return result


def geographic2geocentric(lat):
    # geographic to geocentric
    # https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    # https://instaseis.net/_modules/instaseis/helpers.html
    f = 1.0 / 298.257223563
    e2 = 2 * f - f**2
    lat = np.asarray(lat)
    singulars = (
        (np.abs(lat) < 1e-6) | (np.abs(lat - 90) < 1e-6) | (np.abs(lat + 90) < 1e-6)
    )

    lat_geo = np.rad2deg(np.arctan((1 - e2) * np.tan(np.deg2rad(lat))))
    lat_geo = np.where(singulars, lat, lat_geo)

    # return scalar if input was scalar
    if lat.shape == ():
        return float(lat_geo)
    return lat_geo


def resample_sliprate(slip_rate, dt, dt_new, nsamp):
    """
    For convolution, the sliprate is needed at the sampling of the fields
    in the database. This function resamples the sliprate using linear
    interpolation.

    :param dt: desired sampling
    :param nsamp: desired number of samples
    """
    t_new = np.linspace(0, nsamp * dt_new, nsamp, endpoint=False)
    t_old = np.linspace(0, dt * len(slip_rate), len(slip_rate), endpoint=False)
    return np.interp(t_new, t_old, slip_rate)


def transform_to_spherical(xyz, proj_string, attrs):
    if attrs["coordinates_convention"] == "projected":
        assert (
            proj_string is not None
        ), "coordinates are projected but projection is not defined"
        # Define transformer from projected CRS to geographic
        # (lon, lat, ellipsoidal height)
        transformer = pyproj.Transformer.from_proj(
            pyproj.Proj(proj_string),
            pyproj.Proj(proj="latlong", ellps="sphere", datum="WGS84"),
            always_xy=True,
        )

        lon, lat, _ = transformer.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2])
        xyz[:, 0] = lon
        xyz[:, 1] = lat
        # Depth (xyz[:, 2]) is left unchanged intentionally
    else:
        if attrs["coordinates_convention"] == "geographic":
            xyz[:, 1] = geographic2geocentric(xyz[:, 1])
        elif attrs["coordinates_convention"] == "geocentric":
            print("coordinates already in geocentric")
        else:
            raise ValueError(
                ("coordinates are projected", attrs["coordinates_convention"])
            )
    return xyz


def load_greens_from_hdf5(
    hdf5_file, station_code, fault_tag, segment, rake, components
):
    key = (
        f"{station_code}/"
        f"fault_tag_{fault_tag:03d}/"
        f"segment_{segment[0]:03d}_{segment[1]:03d}/"
        f"rake_{rake}"
    )
    if not os.path.exists(hdf5_file):
        return None

    with h5py.File(hdf5_file, "r") as h5f:
        if key not in h5f:
            return None

        grp = h5f[key]
        traces = []
        for comp in components:
            data = grp[comp][()]
            stats_json = grp.attrs.get(f"{comp}_stats", "{}")
            stats = json.loads(stats_json)
            trace = Trace(data=data)
            trace.stats.update(stats)
            traces.append(trace)
        return Stream(traces)


def save_greens_to_hdf5(
    hdf5_file, station_code, fault_tag, segment, rake, traces, components
):
    key_prefix = (
        f"{station_code}/"
        f"fault_tag_{fault_tag:03d}/"
        f"segment_{segment[0]:03d}_{segment[1]:03d}/"
        f"rake_{rake}"
    )
    with h5py.File(hdf5_file, "a") as h5f:
        grp = h5f.require_group(key_prefix)

        for trace, comp in zip(traces, components):
            # Store data
            if comp in grp:
                del grp[comp]
            grp.create_dataset(comp, data=trace.data)

            # Store stats as JSON string
            stats_str = json.dumps(dict(trace.stats), default=str)
            grp.attrs[f"{comp}_stats"] = stats_str


def check_trace_metadata(trace, src_xyz, lon, lat, epsilon=0.01):
    checks = {
        "source_longitude": (trace.stats.source_longitude, src_xyz[0]),
        "source_latitude": (trace.stats.source_latitude, src_xyz[1]),
        "source_depth_m": (trace.stats.source_depth_m, -src_xyz[2]),
        "station_longitude": (trace.stats.station_longitude, lon),
        "station_latitude": (trace.stats.station_latitude, lat),
    }

    for key, (actual, expected) in checks.items():
        if abs(actual - expected) > epsilon:
            raise ValueError(
                f"Mismatch in {key}: actual={actual:.6f}, expected={expected:.6f}, "
                f"diff={abs(actual - expected):.6f} > epsilon={epsilon}"
            )


def generate_synthetics_instaseis_green(
    db,
    filename,
    t1,
    myproj,
    kind_vd,
    components,
    path_observations,
    station_coords,
    progress_bar,
):
    # read HDF5 and create Finite Source for instaseis
    with h5py.File(filename) as h5f:
        normalized_moment_rates = h5f["normalized_moment_rates"][:, :]
        nsource, ndt = normalized_moment_rates.shape
        xyz = h5f["xyz"][:, :]
        moment_tensors = h5f["moment_tensors"][:, :]
        dt = h5f["dt"][0]
        fault_tags = h5f["fault_tags"][:]
        segment_indices = h5f["segment_indices"][:, :]
        if nsource > 0:
            print(
                f"sources coordinates in {filename}: {xyz[0, :]} ..., {nsource} sources"
            )
        else:
            print(f"{filename} has no sources")
        xyz = transform_to_spherical(xyz, myproj, h5f.attrs)

        # load arguments used for creating the hdf5 file
        json_str = h5f["args_json"][()]
        args_dict = json.loads(json_str)
        dh = int(args_dict["DH"])
        NZ = args_dict["NZ"]
        use_geometric_center = args_dict["use_geometric_center"]
        slip_threshold = args_dict["slip_threshold"]
        assert (
            use_geometric_center
        ), "for using the green function mode, use_geometric_center should be True"
        assert (
            slip_threshold < -1.0
        ), "for using the green function mode, slip_threshold should be small enough"

    db_name = f"{db.info.velocity_model}_{db.info.period}s"
    hdf5_file = f"{path_observations}/greens_{db_name}_dh{dh}_nz{NZ}.h5"
    synthetics = Stream()

    for station_code in station_coords:
        lon, lat = station_coords[station_code]
        network, station = station_code.split(".")

        # create synthetic data with instaseis
        receiver = instaseis.Receiver(
            latitude=geographic2geocentric(lat),
            longitude=lon,
            network=network,
            station=station,
        )
        # print(f"generating instaseis synthetics for station {station}")

        for isrc in range(nsource):
            fault_tag = fault_tags[isrc]
            segment = segment_indices[isrc]
            nodalplane = mt2plane(MomentTensor(moment_tensors[isrc, :], 1))
            M0 = compute_seismic_moment(moment_tensors[isrc, :])
            resampled_moment_rate = resample_sliprate(
                normalized_moment_rates[isrc],
                dt,
                db.info.dt,
                int(ndt * dt / db.info.dt),
            )
            lst = []
            for rake in [0, 90]:
                st0 = load_greens_from_hdf5(
                    hdf5_file, station_code, fault_tag, segment, rake, components
                )
                if st0 is not None:
                    for i in range(len(components)):
                        check_trace_metadata(st0[i], xyz[isrc], lon, lat)

                if not st0:
                    source = instaseis.Source.from_strike_dip_rake(
                        latitude=xyz[isrc, 1],
                        longitude=xyz[isrc, 0],
                        depth_in_m=-xyz[isrc, 2],
                        strike=nodalplane.strike,
                        dip=nodalplane.dip,
                        rake=rake,
                        M0=1.0,
                        origin_time=t1,
                        dt=dt,
                    )
                    sliprate = np.zeros_like(db.info.sliprate)
                    sliprate[1] = 1.0 / db.info.dt
                    st0 = db.get_seismograms(
                        source=source,
                        receiver=receiver,
                        kind=kind_vd,
                        components=components,
                        # sliprate = sliprate
                    )
                    for i in range(len(components)):
                        st0[i].stats.starttime = t1
                        st0[i].stats.source_longitude = xyz[isrc, 0]
                        st0[i].stats.source_latitude = xyz[isrc, 1]
                        st0[i].stats.source_depth_m = -xyz[isrc, 2]
                        st0[i].stats.station_longitude = lon
                        st0[i].stats.station_latitude = lat

                    save_greens_to_hdf5(
                        hdf5_file=hdf5_file,
                        station_code=station_code,
                        fault_tag=fault_tag,
                        segment=segment,
                        rake=rake,
                        components=components,
                        traces=st0,
                    )
                    progress_bar.increment(1)

                lst.append(st0)
            if isrc == 0:
                st = lst[0].copy()
                for i in range(len(components)):
                    st[i].data *= 0
            rake_rad = np.radians(nodalplane.rake)
            for i in range(len(components)):
                G_rake0 = fft_reconvolve_stf(
                    db, lst[0][i].data, np.gradient(resampled_moment_rate)
                )
                G_rake90 = fft_reconvolve_stf(
                    db, lst[1][i].data, np.gradient(resampled_moment_rate)
                )
                st[i].data += M0 * (
                    np.cos(rake_rad) * G_rake0 + np.sin(rake_rad) * G_rake90
                )
        synthetics += st
    return synthetics


def create_finite_source_from_h5(db, filename, t1, myproj):
    # read HDF5 and create Finite Source for instaseis
    with h5py.File(filename) as h5f:
        normalized_moment_rates = h5f["normalized_moment_rates"][:, :]
        nsource, ndt = normalized_moment_rates.shape
        xyz = h5f["xyz"][:, :]
        moment_tensors = h5f["moment_tensors"][:, :]
        dt = h5f["dt"][0]
        print(f"sources coordinates in {filename}", xyz)
        xyz = transform_to_spherical(xyz, myproj, h5f.attrs)

    lps = []
    for isrc in range(nsource):
        source = instaseis.Source(
            latitude=xyz[isrc, 1],
            longitude=xyz[isrc, 0],
            depth_in_m=-xyz[isrc, 2],
            m_rr=moment_tensors[isrc, 0],
            m_tt=moment_tensors[isrc, 1],
            m_pp=moment_tensors[isrc, 2],
            m_rt=moment_tensors[isrc, 3],
            m_rp=moment_tensors[isrc, 4],
            m_tp=moment_tensors[isrc, 5],
            origin_time=t1,
            sliprate=normalized_moment_rates[isrc, :],
            dt=dt,
        )
        source.resample_sliprate(db.info.dt, int(ndt * dt / db.info.dt))
        lps.append(source)
    sources = instaseis.source.FiniteSource(pointsources=lps)
    return sources


def create_finite_source_from_usgs(db, fname, M0_percentile_threshold=0.02):
    sources = instaseis.FiniteSource.from_usgs_param_file(fname, dt=db.info.dt)
    M0 = []
    npts = len(sources.pointsources)
    print(f"usgs file read: {fname}: has {npts} point sources")
    for source in sources.pointsources:
        M0.append(source.M0)
    maxM0 = max(M0)
    filtered_sources = []
    for i, source in enumerate(sources.pointsources):
        if M0[i] > M0_percentile_threshold * maxM0:
            filtered_sources.append(source)
    print(f"filtering low M0 point sources: remaining {len(filtered_sources)}")
    sources = instaseis.source.FiniteSource(pointsources=filtered_sources)
    return sources


class ProgressBarGreen(tqdm):
    """A customized progress bar for tracking the status of synthetics generation."""

    def __init__(self, total):
        super().__init__(total=total)

    def increment(self, number=1):
        self.update(number)


class ProgressBar(tqdm):
    """A customized progress bar for tracking the status of synthetics generation."""

    def __init__(self, total):
        self.total0 = total
        self.current0 = 0
        super().__init__()

    def increment(self, number):
        self.current0 += number

    def update_to(self, current, total):
        self.total = self.total0
        self.update(self.current0 + current - self.n)


def get_unique_segments_per_tag(source_files):
    # Step 1: Collect all available fault tags
    with h5py.File(source_files[0], "r") as h5f:
        fault_tags = h5f["fault_tags"][:]
        available_fault_tags = set(fault_tags)

    print("Available fault tags:", available_fault_tags)

    # Step 2: Dictionary to hold unique segment indices per fault tag
    segments_per_tag = defaultdict(set)

    # Step 3: Loop through files and gather segment indices
    for fname in source_files:
        with h5py.File(fname, "r") as h5f:
            fault_tags = h5f["fault_tags"][:]
            segment_indices = h5f["segment_indices"][:]  # Assumes shape (N, M)

            for k in available_fault_tags:
                ids = np.where(fault_tags == k)[0]  # 1D array of indices
                segs = segment_indices[ids]  # shape (len(ids), M)

                # Convert rows to tuples for hashability and store in set
                for row in segs:
                    segments_per_tag[k].add(tuple(row))

    # Convert sets back to sorted lists of arrays (optional)
    unique_segments_per_tag = {
        k: np.array(sorted(list(v))) for k, v in segments_per_tag.items()
    }
    all_unique_segments = set()

    for segments in unique_segments_per_tag.values():
        for row in segments:
            all_unique_segments.add(tuple(row))  # Convert row to tuple for set

    total_unique_segments = len(all_unique_segments)
    print("Total number of unique segment_indices:", total_unique_segments)

    return total_unique_segments


def generate_synthetics_instaseis_green_function_mode(
    db_name,
    source_files,
    station_coords,
    t1,
    kind_vd,
    components,
    path_observations,
    projection,
):
    db = instaseis.open_db(db_name)

    lst = []
    if len(source_files) == 0:
        return lst

    total_unique_segments = get_unique_segments_per_tag(source_files)
    with ProgressBarGreen(
        2 * total_unique_segments * len(station_coords)
    ) as progress_bar:
        for iModel, fname in enumerate(source_files):
            synth = generate_synthetics_instaseis_green(
                db,
                fname,
                t1,
                projection,
                kind_vd,
                components,
                path_observations,
                station_coords,
                progress_bar,
            )
            lst.append(synth)

    return lst


def generate_synthetics_instaseis_classical_mode(
    db_name,
    source_files,
    station_coords,
    t1,
    kind_vd,
    components,
    path_observations,
    projection,
):
    db = instaseis.open_db(db_name)
    list_finite_sources = []
    for fname in source_files:
        prefix, ex = os.path.splitext(fname)
        if ex == ".param":
            finite_source = create_finite_source_from_usgs(
                db,
                fname,
                M0_percentile_threshold=0.025,
            )
        else:
            finite_source = create_finite_source_from_h5(db, fname, t1, projection)
        list_finite_sources.append(finite_source)

    n_point_sources = [len(sources.pointsources) for sources in list_finite_sources]
    n_total_point_sources = sum(n_point_sources)
    nstations = len(station_coords)
    lst = []
    for iModel, sources in enumerate(list_finite_sources):
        lst.append(Stream())

    with ProgressBar(nstations * n_total_point_sources) as progress_bar:
        for station_code in station_coords:
            lon, lat = station_coords[station_code]
            network, station = station_code.split(".")

            # create synthetic data with instaseis
            receiver = instaseis.Receiver(
                latitude=geographic2geocentric(lat),
                longitude=lon,
                network=network,
                station=station,
            )
            print(f"generating instaseis synthetics for station {station}")

            for iModel, sources in enumerate(list_finite_sources):
                prefix, _ = os.path.splitext(os.path.basename(source_files[iModel]))
                c_time = os.path.getctime(source_files[iModel])
                fname = (
                    f"{path_observations}/{prefix}_{c_time}_{station}_{kind_vd}_"
                    f"{t1.format_iris_web_service()}.mseed"
                )
                if os.path.isfile(fname):
                    print(f"reading the data from {fname}")
                    st0 = read(fname)
                    progress_bar.update_to(n_point_sources[iModel], 0)
                else:
                    print(
                        f"generating synthetics at {network}.{station} for"
                        f" {prefix} ({len(sources)} point sources)"
                    )
                    st0 = db.get_seismograms_finite_source(
                        sources=sources,
                        receiver=receiver,
                        kind=kind_vd,
                        components=components,
                        progress_callback=progress_bar.update_to,
                    )
                for i in range(len(components)):
                    st0[i].stats.starttime = t1
                st0.write(fname, format="MSEED")
                lst[iModel] += st0
                progress_bar.increment(n_point_sources[iModel])
            print("done")
    return lst


def generate_synthetics_instaseis(
    db_name,
    source_files,
    station_coords,
    t1,
    kind_vd,
    components,
    path_observations,
    projection,
    modes=["classical"],
):
    lst = []
    if "classical" in modes:
        lst1 = generate_synthetics_instaseis_classical_mode(
            db_name,
            source_files,
            station_coords,
            t1,
            kind_vd,
            components,
            path_observations,
            projection,
        )
        lst += lst1

    if "green_functions" in modes:
        lst2 = generate_synthetics_instaseis_green_function_mode(
            db_name,
            source_files,
            station_coords,
            t1,
            kind_vd,
            components,
            path_observations,
            projection,
        )
        lst += lst2
    if ("green_functions" not in modes) and ("classical" not in modes):
        raise ValueError(f"unknown instaseis modes {modes}")
    return lst
