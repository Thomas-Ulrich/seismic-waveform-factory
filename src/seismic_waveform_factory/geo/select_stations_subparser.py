import argparse


def add_parser(subparsers):
    parser = subparsers.add_parser(
        "select-stations",
        help="Select stations ensuring optimal coverage from configuration file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional arguments
    parser.add_argument("config_file", help="yaml config file describing the event.")
    parser.add_argument(
        "number_stations", type=int, help="Total number of stations to select."
    )
    parser.add_argument(
        "closest_stations", type=int, help="Number of closest stations to select."
    )

    # Optional arguments
    parser.add_argument(
        "--distance_range",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="Distance range (degrees) from which to select stations.",
    )
    parser.add_argument(
        "--channel",
        default="*",
        type=str,
        help='Filter channels to be retrieved (default "*" = all channels).',
    )
    parser.add_argument(
        "--store_format",
        choices=["sac", "mseed"],
        default="mseed",
        type=str,
        help="Storage format for waveform data.",
    )
    parser.add_argument(
        "--azimuthal",
        action="store_true",
        help="Select stations based on back-azimuth instead of distance.",
    )
    parser.add_argument(
        "--station_kind",
        choices=["auto", "regional", "global"],
        default="auto",
        type=str,
        help="kind of stations to select.",
    )

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.geo.select_stations import main

        main(args)

    parser.set_defaults(func=run)
