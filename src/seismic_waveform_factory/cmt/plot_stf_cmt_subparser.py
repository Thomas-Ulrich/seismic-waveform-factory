import argparse


def add_parser(subparsers):
    parser = subparsers.add_parser(
        "plot-stf-cmt",
        help="Plot source time functions from a multi-CMT file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("filename", help="Moment Tensor (h5)")
    parser.add_argument(
        "--idSTF", nargs="+", help="list of the STF to visualize (1st = 0); ", type=int
    )
    parser.add_argument(
        "--normalized",
        dest="normalized",
        action="store_true",
        help="show normalized stf (integral=1)",
    )
    parser.add_argument(
        "--cumulative",
        dest="cumulative",
        action="store_true",
        help="show a cumulative plot of all STF",
    )
    parser.add_argument(
        "--MRF", nargs=1, metavar=("MomentRateFile"), help="save MomentRate to file"
    )
    parser.add_argument(
        "--time_range",
        nargs=2,
        metavar=("minT", "maxT"),
        help="time range for time slices",
        type=float,
    )
    parser.add_argument(
        "--extension", nargs=1, default=(["svg"]), help="extension output file"
    )

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.cmt.plot_stf_cmt import main

        main(args)

    parser.set_defaults(func=run)
