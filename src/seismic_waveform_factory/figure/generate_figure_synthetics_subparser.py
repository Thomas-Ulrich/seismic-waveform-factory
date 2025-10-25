def add_parser(subparsers):
    parser = subparsers.add_parser(
        "plot-waveforms", help="Generate waveform figures from configuration file"
    )
    parser.add_argument("config_file", help="configuration file")

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.figure.generate_figure_syntetics import main

        main(args)

    parser.set_defaults(func=run)
