def add_parser(subparsers):
    parser = subparsers.add_parser(
        "gen-map-waveforms", help="Generate map with waveforms from configuration file"
    )
    parser.add_argument("config_file", help="config file describing event and stations")
    parser.add_argument(
        "plot_panel", help="panel from which to extract waveforms", type=int
    )

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.figure.generate_map_with_waveforms import main

        main(args)

    parser.set_defaults(func=run)
