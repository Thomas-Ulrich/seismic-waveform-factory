import argparse


def add_parser(subparsers):
    parser = subparsers.add_parser(
        "plot-waveforms",
        help="Generate waveform figures from configuration file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("config_file", help="configuration file")

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.figure.generate_figure_synthetics import main

        main(args)

    parser.set_defaults(func=run)
