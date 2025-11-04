import argparse


def add_parser(subparsers):
    parser = subparsers.add_parser(
        "plot-station-map",
        help="Plot station map from configuration file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("config_file", help="config file describing event and stations")
    parser.add_argument("--plot_all_station_file", action="store_true")

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.figure.plot_station_map import main

        main(args)

    parser.set_defaults(func=run)
