import argparse


def add_parser(subparsers_top):
    parser = subparsers_top.add_parser(
        "compute-multi-cmt",
        help="Compute multi-CMT representation from SeisSol fault output file",
    )

    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True

    spatial = subparsers.add_parser(
        "spatial",
        help="Divide fault along horizontal vector Vh and Uz",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
            action="store_true",
            help=("compute potency instead of seismic moment (basically use G=1)"),
        )

        sp.add_argument(
            "--STFfromSR",
            metavar=("xdmf SR File"),
            help=(
                "Use SR to compute Source Time Function "
                "(high sampling rate required)"
            ),
        )

        sp.add_argument(
            "--proj",
            metavar=("projname"),
            help=(
                "Transform to the coordinate reference system WGS 84 (lat, lon) "
                "from projname"
            ),
        )

        sp.add_argument(
            "--DH",
            default=20,
            type=float,
            help="Max horizontal distance between point sources, in km",
        )

        sp.add_argument(
            "--NZ",
            metavar=("nz"),
            default=2,
            type=int,
            help="Number of point sources along z",
        )

        sp.add_argument(
            "--slip_threshold",
            metavar=("slip_threshold"),
            default=0.1,
            type=float,
            help=("Slip threshold for excluding low slip areas when slicing the fault"),
        )

        sp.add_argument(
            "--refVector",
            nargs=3,
            metavar=("rvx", "rvy", "rvz"),
            default=([-1e-5, 0, -1]),
            type=float,
            help="Reference vector as defined in SeisSol, for computing strike and dip",
        )

        sp.add_argument(
            "--ndt",
            metavar=("ndt"),
            type=int,
            help="Use a subset of time frames",
        )

        sp.add_argument(
            "--invertSld",
            action="store_true",
            help=(
                "Invert Sld (if the normal is consistently wrongly defined "
                "in SeisSol output)"
            ),
        )

        sp.add_argument(
            "--use_geometric_center",
            action="store_true",
            help=(
                "Use the geometric center (area-weighted centroid) of the selected "
                " face instead of the moment-weighted center to define the point source"
                " location. This makes the location independent of the moment "
                "distribution, which is useful when reusing Green's functions across "
                "different rupture models."
            ),
        )

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.cmt.compute_multi_cmt import main

        main(args)

    parser.set_defaults(func=run)
