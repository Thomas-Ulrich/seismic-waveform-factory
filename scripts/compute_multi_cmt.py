#!/usr/bin/env python3
import numpy as np
import cmt
import argparse
from faultoutput import FaultOutput
import os

parser = argparse.ArgumentParser(
    description=(
        "Compute multiple point source representation from SeisSol fault output"
    )
)
subparsers = parser.add_subparsers(dest="command")
subparsers.required = True

spatial = subparsers.add_parser(
    "spatial",
    help="Divide fault along horizontal vector Vh and Uz",
)

spatial.add_argument(
    "--vH",
    nargs=2,
    metavar=("vhx", "vhy"),
    help="Vector defining the slicing direction in the horizontal plane",
    type=float,
)

temporal = subparsers.add_parser(
    "temporal",
    help="Divide fault with time and Uz",
)

temporal.add_argument(
    "--time_range",
    nargs=2,
    metavar=("minT", "maxT"),
    help="Time range for time slices",
    type=float,
)

temporal.add_argument(
    "--time_sub_events",
    nargs="+",
    type=float,
    help=(
        "When temporal, define the exact times (relative to minT, see time_range) "
        "separating events (NH is then ignored)"
    ),
)

for sp in [spatial, temporal]:
    sp.add_argument(
        "filename",
        help="Fault output filename (XDMF)",
    )

    sp.add_argument(
        "mu",
        help=(
            "a float value, a 2 columns ASCII (z, mu) text file "
            "or a yaml file to be read with easi"
        ),
    )

    sp.add_argument(
        "--potency",
        dest="potency",
        action="store_true",
        help=("compute potency instead of seismic moment (basically use G=1)"),
    )

    sp.add_argument(
        "--STFfromSR",
        nargs=1,
        metavar=("xdmf SR File"),
        help=(
            "Use SR to compute Source Time Function " "(high sampling rate required)"
        ),
    )

    sp.add_argument(
        "--proj",
        nargs=1,
        metavar=("projname"),
        help=(
            "Transform to the coordinate reference system WGS 84 (lat, lon) "
            "from projname"
        ),
    )

    sp.add_argument(
        "--DH",
        nargs=1,
        default=([20]),
        type=float,
        help="Max horizontal distance between point sources, in km",
    )

    sp.add_argument(
        "--NZ",
        nargs=1,
        metavar=("nz"),
        default=([2]),
        type=int,
        help="Number of point sources along z",
    )

    sp.add_argument(
        "--slip_threshold",
        nargs=1,
        metavar=("slip_threshold"),
        default=([0.1]),
        type=float,
        help=(
            "Slip threshold used for excluding low slip areas when slicing the fault"
        ),
    )

    sp.add_argument(
        "--refVector",
        nargs=3,
        metavar=("rvx", "rvy", "rvz"),
        default=([-1e-5, 0, -1]),
        type=float,
        help="Reference vector as defined in SeisSol (for computing strike and dip)",
    )

    sp.add_argument(
        "--ndt",
        nargs=1,
        metavar=("ndt"),
        type=int,
        help="Use a subset of time frames",
    )

    sp.add_argument(
        "--invertSld",
        dest="invertSld",
        action="store_true",
        help=(
            "Invert Sld (if the normal is consistently wrongly defined "
            "in SeisSol output)"
        ),
    )

args = parser.parse_args()

if args.proj:
    from pyproj import Transformer

    # epsg:4326 is the coordinate reference system WGS 84 (lat, lon)
    transformer = Transformer.from_crs(args.proj[0], "epsg:4326", always_xy=True)

fo = FaultOutput(args.filename, args.ndt)
fo.read_final_slip()
fo.compute_strike_dip(args.refVector)
fo.compute_rake(args.invertSld)
fo.compute_barycenter_coords()
fo.compute_face_area()

if args.potency:
    fo.G = 1.0
else:
    fo.evaluate_G(args.mu)

fo.Garea = fo.G * fo.face_area

if args.STFfromSR:
    fo_SR = FaultOutput(args.STFfromSR[0])
    fo_SR.compute_barycenter_coords()
    fo.compute_face_moment_rate_from_slip_rate(fo_SR)
else:
    fo.compute_face_moment_rate_from_ASl()

fo.compute_face_moment_tensor_NED()

if args.command == "temporal":
    fo.fault_tags[:] = 3
    fo.unique_fault_tags = [3]

M0_eq = 0
nsrc_eq = 0
point_sources = {}

for fault_tag in fo.unique_fault_tags:
    selected = np.where(fo.fault_tags == fault_tag)[0]
    if args.command == "temporal":
        if args.time_range:
            t0 = args.time_range[0]
            tmax = args.time_range[1]
        else:
            t0 = 0.0
            tmax = fo.dt * fo.ndt
        if args.time_sub_events:
            t_abs_slices = [t0 + val for val in args.time_sub_events]
            hslices = [t0, *t_abs_slices, tmax]
            sprint = f"user defined time for sub_events: {hslices}"
        else:
            hslices = np.linspace(t0, tmax, (args.DH[0] + 1))
            sprint = f"{len(hslices) - 1} slices along rupture time ({t0} - {tmax}s)"
        # hslices[-1] = 1e10
        h_or_RT = fo.get_rupture_time()
    else:
        if not args.vH:
            from sklearn.decomposition import PCA

            print("args.VH not set... inferring from PCA")
            # Perform PCA to get principal axes
            pca = PCA(n_components=2)
            points = pca.fit_transform(fo.xyzc[selected, 0:2])
            ua, ub = pca.components_
            print(ua, ub)
            args.vH = [ua[0], ua[1]]

        h_or_RT = args.vH[0] * fo.xyzc[selected, 0] + args.vH[1] * fo.xyzc[selected, 1]
        hslices = cmt.compute_slices_array_enforcing_dx(
            h_or_RT, fo.slip[selected], args.DH[0], args.slip_threshold[0]
        )
        sprint = f"{len(hslices) - 1} slices along horizontal direction :\
         ({args.vH[0]},{args.vH[1]})"

    print(f"{sprint}, {args.NZ[0]} along z")
    zcenters0 = fo.xyzc[selected, 2]
    zslices = cmt.compute_slices_array(
        zcenters0, fo.slip[selected], args.NZ[0], args.slip_threshold[0]
    )

    nsources = (len(hslices) - 1) * args.NZ[0]
    aNormMRF = np.zeros((nsources, fo.ndt))
    aMomentTensor = np.zeros((nsources, 6))
    axyz = np.zeros((nsources, 3))

    M0_segment = 0
    isrc_segment = 0
    for i in range(len(hslices) - 1):
        h_or_RT1, h_or_RT2 = hslices[i], hslices[i + 1]
        for j in range(args.NZ[0]):
            z1, z2 = zslices[j], zslices[j + 1]

            idxys = np.where((h_or_RT >= h_or_RT1) & (h_or_RT < h_or_RT2))
            idzs = np.where((zcenters0 >= z1) & (zcenters0 < z2))
            ids = np.intersect1d(idxys[0], idzs[0])

            (
                M0,
                NormMRF,
                MomentTensor,
                xyz,
            ) = fo.compute_equivalent_point_source_subfault(selected[ids])
            if not M0:
                continue

            if args.proj:
                # project back to geocentric coordinates
                xyz[0], xyz[1] = transformer.transform(xyz[0], xyz[1])

            M0_segment = M0_segment + M0
            (
                aNormMRF[isrc_segment, :],
                aMomentTensor[isrc_segment, :],
                axyz[isrc_segment, :],
            ) = (
                NormMRF,
                MomentTensor,
                xyz,
            )
            isrc_segment = isrc_segment + 1

    M0_eq += M0_segment
    nsrc_eq += isrc_segment
    if args.potency:
        print(
            f"Potency(segment_{fault_tag})={M0_segment:.2e} m3"
            f", {isrc_segment} sources"
        )
    else:
        Mw = 2.0 / 3.0 * np.log10(M0_segment) - 6.07
        print(
            f"Mw(segment_{fault_tag})={Mw:.2f} ({M0_segment:.2e} Nm)"
            f", {isrc_segment} sources"
        )

    aMomentTensor = cmt.NED2RTP(aMomentTensor[0:isrc_segment, :])
    aNormMRF = aNormMRF[0:isrc_segment, :]
    axyz = axyz[0:isrc_segment, :]

    point_sources[fault_tag] = {
        "moment_tensors": aMomentTensor,
        "moment_rate_functions": aNormMRF,
        "locations": axyz,
    }

prefix = os.path.basename(args.filename.split("-fault")[0])
if args.command == "temporal":
    fname = f"PointSourceFile_{prefix}_nt{args.DH[0]}_nz{args.NZ[0]}.h5"
else:
    fname = f"PointSourceFile_{prefix}_dx{args.DH[0]}_nz{args.NZ[0]}.h5"
cmt.write_point_source_file(fname, point_sources, fo.dt, args.proj, args.potency)

if args.potency:
    print(f"Potency(earthquake)= {M0_eq:.2e} m3, {nsrc_eq} sources")
else:
    Mw = 2.0 / 3.0 * np.log10(M0_eq) - 6.07
    print(f"Mw(earthquake)={Mw:.2f} ({M0_eq:.2e} Nm), {nsrc_eq} sources")
