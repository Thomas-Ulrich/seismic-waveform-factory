import instaseis
from tqdm import tqdm
import h5py
import os
from obspy import read, Stream


def geographic2geocentric(lat):
    import numpy as np

    # geographic to geocentric
    # https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    f = 1.0 / 298.257223563
    e2 = 2 * f - f**2
    return np.rad2deg(np.arctan((1 - e2) * np.tan(np.deg2rad(lat))))


def create_finite_source_from_h5(db, filename, t1, myproj):
    # read HDF5 and create Finite Source for instaseis
    with h5py.File(filename) as h5f:
        NormalizedMomentRate = h5f["NormalizedMomentRate"][:, :]
        nsource, ndt = NormalizedMomentRate.shape
        xyz = h5f["xyz"][:, :]
        MomentTensor = h5f["MomentTensor"][:, :]
        dt = h5f["dt"][0]
        print(f"sources coordinates in {filename}", xyz)
        if myproj:
            print("projecting back to geocentric")
            import pyproj

            myproj = pyproj.Proj(myproj)
            lla = pyproj.Proj(proj="latlong", ellps="sphere", datum="WGS84")
            txyz = pyproj.transform(
                myproj, lla, xyz[:, 0], xyz[:, 1], xyz[:, 2], radians=False
            )
            xyz[:, 0] = txyz[0]
            xyz[:, 1] = txyz[1]
            # we use the same depth (else it is modified by the ellipsoid to sphere)
            # xyz[:,2] = txyz[2]
            print(xyz)
        else:
            if h5f.attrs["CoordinatesConvention"] == b"geographic":
                xyz[:, 1] = geographic2geocentric(xyz[:, 1])
            elif h5f.attrs["CoordinatesConvention"] == b"geocentric":
                print("coordinates already in geocentric")
            else:
                print("coordinates are projected", h5f.attrs["CoordinatesConvention"])
                exit()
    lps = []
    for isrc in range(nsource):
        source = instaseis.Source(
            latitude=xyz[isrc, 1],
            longitude=xyz[isrc, 0],
            depth_in_m=-xyz[isrc, 2],
            m_rr=MomentTensor[isrc, 0],
            m_tt=MomentTensor[isrc, 1],
            m_pp=MomentTensor[isrc, 2],
            m_rt=MomentTensor[isrc, 3],
            m_rp=MomentTensor[isrc, 4],
            m_tp=MomentTensor[isrc, 5],
            origin_time=t1,
            sliprate=NormalizedMomentRate[isrc, :],
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


def generate_synthetics_instaseis(
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
