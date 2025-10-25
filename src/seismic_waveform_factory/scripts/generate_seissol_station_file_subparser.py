def add_parser(subparsers):
    parser = subparsers.add_parser(
        "gen-seissol-sta", help="Generate SeisSol station file from configuration file"
    )
    parser.add_argument("config_file", help="configuration file")

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.scripts.generate_seissol_station_file import main

        main(args)

    parser.set_defaults(func=run)
