#!/usr/bin/env python3
import json
import os
import copy

import numpy as np

from seismic_waveform_factory.fault.fault_output import FaultOutput
from seismic_waveform_factory.utils import cmt


def main(args):
    if args.proj:
        from pyproj import Transformer

        # epsg:4326 is the coordinate reference system WGS 84 (lat, lon)
        transformer = Transformer.from_crs(args.proj, "epsg:4326", always_xy=True)

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
        fo_SR = FaultOutput(args.STFfromSR)
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
                hslices = np.linspace(t0, tmax, (args.DH + 1))
                sprint = (
                    f"{len(hslices) - 1} slices along rupture time ({t0} - {tmax}s)"
                )
            # hslices[-1] = 1e10
            h_or_RT = fo.get_rupture_time()
        else:
            if not args.vH:
                from sklearn.decomposition import PCA

                print("args.VH not set... inferring from PCA")
                # Perform PCA to get principal axes
                pca = PCA(n_components=2)
                pca.fit_transform(fo.xyzc[selected, 0:2])
                ua, ub = pca.components_
                print(ua, ub)
                args.vH = [ua[0], ua[1]]

            h_or_RT = (
                args.vH[0] * fo.xyzc[selected, 0] + args.vH[1] * fo.xyzc[selected, 1]
            )
            hslices = cmt.compute_slices_array_enforcing_dx(
                h_or_RT, fo.slip[selected], args.DH, args.slip_threshold
            )
            sprint = f"{len(hslices) - 1} slices along horizontal direction :\
             ({args.vH[0]},{args.vH[1]})"

        print(f"{sprint}, {args.NZ} along z")
        zcenters0 = fo.xyzc[selected, 2]
        zslices = cmt.compute_slices_array(
            zcenters0, fo.slip[selected], args.NZ, args.slip_threshold
        )

        nsources = (len(hslices) - 1) * args.NZ
        aNormMRF = np.zeros((nsources, fo.ndt))
        aMomentTensor = np.zeros((nsources, 6))
        axyz = np.zeros((nsources, 3))
        segment_indices = np.zeros((nsources, 2), dtype=int)

        M0_segment = 0
        isrc_segment = 0
        for i in range(len(hslices) - 1):
            h_or_RT1, h_or_RT2 = hslices[i], hslices[i + 1]
            for j in range(args.NZ):
                z1, z2 = zslices[j], zslices[j + 1]

                idxys = np.where((h_or_RT >= h_or_RT1) & (h_or_RT < h_or_RT2))
                idzs = np.where((zcenters0 >= z1) & (zcenters0 < z2))
                ids = np.intersect1d(idxys[0], idzs[0])

                (
                    M0,
                    NormMRF,
                    MomentTensor,
                    xyz,
                ) = fo.compute_equivalent_point_source_subfault(
                    selected[ids], use_geometric_center=args.use_geometric_center
                )
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
                segment_indices[isrc_segment, :] = i, j
                isrc_segment = isrc_segment + 1

        M0_eq += M0_segment
        nsrc_eq += isrc_segment
        if args.potency:
            print(
                f"Potency(segment_{fault_tag})={M0_segment:.2e} m3"
                f", {isrc_segment} sources"
            )
        else:
            Mw = 2.0 / 3.0 * np.log10(M0_segment) - 6.07 if M0_segment else 0.0
            print(
                f"Mw(segment_{fault_tag})={Mw:.2f} ({M0_segment:.2e} Nm)"
                f", {isrc_segment} sources"
            )

        aMomentTensor = cmt.NED2RTP(aMomentTensor[0:isrc_segment, :])
        aNormMRF = aNormMRF[0:isrc_segment, :]
        axyz = axyz[0:isrc_segment, :]
        segment_indices = segment_indices[0:isrc_segment, :]

        point_sources[fault_tag] = {
            "moment_tensors": aMomentTensor,
            "moment_rate_functions": aNormMRF,
            "locations": axyz,
            "segment_indices": segment_indices,
        }
    # Deep copy to avoid modifying the real one
    args_copy = copy.deepcopy(vars(args))
    # Remove non-serializable fields
    args_copy.pop("func", None)
    json_str = json.dumps(args_copy)

    prefix = os.path.basename(args.filename.split("-fault")[0])
    if args.command == "temporal":
        fname = f"PointSourceFile_{prefix}_nt{args.DH}_nz{args.NZ}.h5"
    else:
        fname = f"PointSourceFile_{prefix}_dx{args.DH}_nz{args.NZ}.h5"
    cmt.write_point_source_file(
        fname, point_sources, fo.dt, args.proj, args.potency, json_str
    )

    if args.potency:
        print(f"Potency(earthquake)= {M0_eq:.2e} m3, {nsrc_eq} sources")
    else:
        Mw = 2.0 / 3.0 * np.log10(M0_eq) - 6.07
        print(f"Mw(earthquake)={Mw:.2f} ({M0_eq:.2e} Nm), {nsrc_eq} sources")
