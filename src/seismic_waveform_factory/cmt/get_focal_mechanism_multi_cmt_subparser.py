def add_parser(subparsers):
    parser = subparsers.add_parser(
        "get-focal-mech", help="Extract focal mechanism from multi-CMT file"
    )
    parser.add_argument("filename", help="point source file (h5)")

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.cmt.get_focal_mechanism_multi_cmt import main

        main(args)

    parser.set_defaults(func=run)
