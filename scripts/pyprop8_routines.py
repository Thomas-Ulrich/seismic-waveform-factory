import sys
import numpy as np
from obspy import UTCDateTime
from obspy import read, Trace
from obspy.imaging.beachball import mt2plane, MomentTensor
import cmt
from pyproj import Transformer
import h5py
from scipy.signal import fftconvolve
import pyprop8 as pp
from pyprop8.utils import stf_trapezoidal, make_moment_tensor, rtf2xyz
from tqdm import tqdm
import os


def create_pyprop8_receivers(list_inventory, transformer):
    station_coords = np.zeros((len(list_inventory), 2))
    for ins, inv in enumerate(list_inventory):
        sta = inv[0][0]
        station_coords[ins, :] = [sta.longitude, sta.latitude]
    if transformer:
        station_coords[:, 0], station_coords[:, 1] = transformer.transform(
            station_coords[:, 0], station_coords[:, 1]
        )
        station_coords[:, :] /= 1000.0
        geometry = "cartesian"
    else:
        geometry = "spherical"
    print(f"receivers ({geometry})", station_coords)
    receivers = pp.ListOfReceivers(
        station_coords[:, 0], station_coords[:, 1], 0, geometry=geometry
    )
    return receivers


def create_pyprop8_source_list_from_h5(filename):
    with h5py.File(filename) as h5f:
        normalizedMomentRate = h5f["NormalizedMomentRate"][:, :]
        nsource, ndt = normalizedMomentRate.shape
        xyz = h5f["xyz"][:, :]
        aMomentTensor = h5f["MomentTensor"][:, :]
        dt = h5f["dt"][0]
        assert h5f.attrs["CoordinatesConvention"] == b"geographic"

    latlon = True
    if not latlon:
        lat0, lon0 = xyz[:, 1].mean(), xyz[:, 0].mean()
        proj = f"+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0={lon0} +lat_0={lat0}"
        transformer = Transformer.from_crs("epsg:4326", proj, always_xy=True)
        xyz[:, 0], xyz[:, 1] = transformer.transform(xyz[:, 0], xyz[:, 1])
        xyz[:, :] /= 1000.0
    else:
        transformer = None
        xyz[:, 2] /= 1000.0

    M = np.zeros((3, 3))
    F = np.zeros((3, 1))
    sources = []
    for isrc in range(nsource):
        M[0, 0] = aMomentTensor[isrc, 2]
        M[1, 1] = aMomentTensor[isrc, 1]
        M[2, 2] = aMomentTensor[isrc, 0]
        M[0, 1] = M[1, 0] = -aMomentTensor[isrc, 5]
        M[0, 2] = M[2, 0] = aMomentTensor[isrc, 4]
        M[1, 2] = M[2, 1] = -aMomentTensor[isrc, 3]
        M[:, :] *= 1e-15

        use_strike_dip_rake = False
        if use_strike_dip_rake:
            M0all = cmt.compute_seismic_moment(aMomentTensor[isrc, :])
            # Mw = 2.0 /    3.0 * np.log10(M0all) - 6.07
            nopl = mt2plane(MomentTensor(aMomentTensor[isrc, :], 0))
            M[:, :] = rtf2xyz(
                make_moment_tensor(
                    nopl.strike, nopl.dip, nopl.rake, M0all * 1e-15, 0, 0
                )
            )

        sources += [
            pp.PointSource(
                xyz[isrc, 0],
                xyz[isrc, 1],
                -xyz[isrc, 2],
                M[:, :],
                F[:],
                0.0,
            )
        ]

    return transformer, sources, dt, normalizedMomentRate


def generate_synthetics(
    dt_stf, stf, nstation, kind_vd, model, source, receivers, fmax, duration
):
    nsrc, ndt = stf.shape
    time_stf = np.arange(ndt) * dt_stf

    dt = 1 / fmax
    nt = int(duration / dt)
    time_resampled = np.arange(0, time_stf[-1], dt)
    synth = np.zeros((nstation * 3, nt))
    try:
        num_proc = os.environ["OMP_NUM_THREADS"]
        print(f"using {num_proc} thread in pyprop8")
    except KeyError:
        print("OMP_NUM_THREADS not set, using 1 thread in pyprop8")
        num_proc = 1
    for isrc in tqdm(range(nsrc)):
        t, green_functions = pp.compute_seismograms(
            model,
            source[isrc],
            receivers,
            nt,
            dt,
            squeeze_outputs=False,
            number_of_processes=8,
        )

        resampled_stf = np.interp(time_resampled, time_stf, stf[isrc, :])

        assert np.abs(dt - t[1]) < 1e-6
        for ista in range(nstation):
            for k in range(3):
                synth[k * nstation + ista] += fftconvolve(
                    dt * green_functions[0, ista, k, :], resampled_stf, "full"
                )[0:nt]

    if kind_vd in ["velocity", "accelaration"]:
        synth = np.gradient(synth, axis=-1) / dt
    if kind_vd == "accelaration":
        synth = np.gradient(synth, axis=-1) / dt

    return dt, synth


def create_synthetic_stream(starttime, dt, list_inventory, synth):
    st_syn = read()
    st_syn.clear()
    xyz = "ENZ"
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


def generate_synthetics_pyprop8(
    source_files,
    list_inventory,
    onset,
    kind_vd,
    fmax,
    duration,
    vel_model_fname,
):
    lst = []

    model_axitra_fmt = np.loadtxt(vel_model_fname, comments="#")
    model_axitra_fmt = model_axitra_fmt[:, 0:4] / 1e3
    if model_axitra_fmt[-1, 0] == 0:
        model_axitra_fmt[-1, 0] = np.inf
    model = pp.LayeredStructureModel(model_axitra_fmt)
    print(model)
    nreceivers = len(list_inventory)

    for fname in source_files:
        transformer, lsources, dt_stf, stf = create_pyprop8_source_list_from_h5(fname)
        receivers = create_pyprop8_receivers(list_inventory, transformer)
        dt, synth = generate_synthetics(
            dt_stf, stf, nreceivers, kind_vd, model, lsources, receivers, fmax, duration
        )
        st_syn = create_synthetic_stream(onset, dt, list_inventory, synth)
        lst.append(st_syn)
    return lst
