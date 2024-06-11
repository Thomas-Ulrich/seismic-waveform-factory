import sys
import numpy as np
from obspy import UTCDateTime
from obspy import read, Trace
from obspy.imaging.beachball import mt2plane, MomentTensor
import cmt
from pyproj import Transformer
from axitra import Axitra, moment
import h5py
from scipy.signal import fftconvolve
from tqdm import tqdm


def create_axitra_station_file(list_inventory):
    stations_array_axitra = np.zeros((len(list_inventory), 4))
    for ins, inv in enumerate(list_inventory):
        sta = inv[0][0]
        stations_array_axitra[ins, :] = [ins + 1, sta.latitude, sta.longitude, 0]
    return stations_array_axitra


def create_axitra_source_from_h5(filename):
    with h5py.File(filename) as h5f:
        normalizedMomentRate = h5f["NormalizedMomentRate"][:, :]
        nsource, ndt = normalizedMomentRate.shape
        xyz = h5f["xyz"][:, :]
        aMomentTensor = h5f["MomentTensor"][:, :]
        dt = h5f["dt"][0]
        sources = np.zeros((nsource, 4))
        sources[:, 0] = np.arange(nsource) + 1
        sources[:, 1] = xyz[:, 1]
        sources[:, 2] = xyz[:, 0]
        sources[:, 3] = -xyz[:, 2]

        assert h5f.attrs["CoordinatesConvention"] == b"geographic"

    hist = np.zeros((nsource, 8))
    delay = 0
    for isrc in range(nsource):
        M0all = cmt.compute_seismic_moment(aMomentTensor[isrc, :])
        # Mw = 2.0 / 3.0 * np.log10(M0all) - 6.07
        nopl = mt2plane(MomentTensor(aMomentTensor[isrc, :], 0))
        hist[isrc, :] = [
            isrc + 1,
            M0all,
            nopl.strike,
            nopl.dip,
            nopl.rake,
            0,
            0,
            delay,
        ]
    return sources, hist, dt, normalizedMomentRate


def compute_dt_nt(ap):
    nt = int(2 ** (1 + np.ceil(np.log2((ap.duration - 0.499999) * ap.fmax))))
    dt = ap.duration / nt
    return dt, nt


def generate_synthetics_from_axitra_green_functions(
    ap, dt_stf, stf, hist, nstation, kind_vd
):
    nsrc, ndt_stf = stf.shape
    time_stf = np.arange(ndt_stf) * dt_stf

    dt_axitra, ndt_axitra = compute_dt_nt(ap)
    time_resampled = np.arange(0, time_stf[-1], dt_axitra)
    synth = np.zeros((nstation * 3, ndt_axitra))
    output_dic = {"acceleration": 3, "velocity": 2, "displacement": 1}

    for isrc in range(nsrc):
        histc = hist.copy()
        # we cant convolve with a set of source time function for each source
        # therefore we convolve one by one, and delay the other sources
        histc[:, -1] = 1e10
        histc[isrc, -1] = 0
        resampled_stf = np.interp(time_resampled, time_stf, stf[isrc, :])
        # 4 seems to decrease the HF noise
        # we convolve manually because I cant
        # make sfunc work properly
        t, sx, sy, sz = moment.conv(
            ap,
            histc,
            source_type=4,
            t0=dt_axitra,
            unit=output_dic[kind_vd],
        )
        assert np.abs(dt_axitra - t[1]) < 1e-6
        gf_list = [sx, sy, sz]
        for ista in range(nstation):
            for k, gf in enumerate(gf_list):
                # the dt_axitra normalizes the dirac
                synth[k * nstation + ista] += fftconvolve(
                    gf[ista, :], resampled_stf, "full"
                )[0:ndt_axitra]
    # 0.01 : cm2m 0.5*dt_axitra: normalization of the Dirac?
    return dt_axitra, synth * dt_axitra * 0.5 * 0.01


def create_synthetic_stream(starttime, dt, list_inventory, synth):
    st_syn = read()
    st_syn.clear()
    xyz = "NEZ"
    nsta = len(list_inventory)
    for ista, sta in enumerate(list_inventory):
        for i in range(0, 3):
            tr = Trace()
            tr.stats.network = sta[0].code
            tr.stats.station = sta[0][0].code
            tr.stats.channel = xyz[i]
            tr.data = synth[i * nsta + ista, :]
            tr.stats.delta = dt
            tr.stats.starttime = starttime
            st_syn.append(tr)
    return st_syn


def generate_synthetics_axitra(
    source_files,
    list_inventory,
    onset,
    kind_vd,
    fmax,
    duration,
    vel_model_fname,
    path_axitra,
):
    lst = []
    nsta = len(list_inventory)

    for fname in tqdm(source_files):
        # because of ap.clean, model and station_array_axitra need to be rewritten
        model = np.loadtxt(vel_model_fname, comments="#")
        stations_array_axitra = create_axitra_station_file(list_inventory)
        sources, hist, dt_stf, stf = create_axitra_source_from_h5(fname)
        latlon = False
        if not latlon:
            lat0, lon0 = sources[:, 1].mean(), sources[:, 2].mean()
            proj = f"+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0={lon0} +lat_0={lat0}"
            transformer = Transformer.from_crs("epsg:4326", proj, always_xy=False)

            # x=north, y=east, z=upward therefore the index swap
            arrays = [sources, stations_array_axitra]
            for array in arrays:
                array[:, 2], array[:, 1] = transformer.transform(
                    array[:, 1], array[:, 2]
                )

        ap = Axitra(
            model,
            stations_array_axitra,
            sources,
            fmax=fmax,
            duration=duration,
            xl=0.0,
            latlon=latlon,
            axpath=path_axitra,
        )
        ap = moment.green(ap)

        dt_axitra, synth = generate_synthetics_from_axitra_green_functions(
            ap, dt_stf, stf, hist, nsta, kind_vd
        )
        st_syn = create_synthetic_stream(onset, dt_axitra, list_inventory, synth)
        lst.append(st_syn)
        ap.clean()
    return lst
