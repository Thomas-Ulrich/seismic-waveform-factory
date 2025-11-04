import argparse


def add_parser(subparsers):
    parser = subparsers.add_parser(
        "gen-legend-box",
        help="Generate legend box from configuration file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("config_file", help="config file describing event and stations")

    def run(args):
        from seismic_waveform_factory.figure.generate_legend_box import main

        main(args)

    parser.set_defaults(func=run)
