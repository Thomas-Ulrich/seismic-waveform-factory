def add_parser(subparsers):
    parser = subparsers.add_parser(
        "plot-map-cmt", help="Plot a map including beachballs from multi-CMT file"
    )

    parser.add_argument("filename", help="point source file (h5)")
    parser.add_argument(
        "--proj",
        nargs=1,
        metavar=("projname"),
        help="transform source coordinates from projection projname to geocentric",
    )
    parser.add_argument(
        "--x0y0proj",
        nargs=2,
        metavar=("x0", "y0"),
        default=[0, 0],
        help="offset of projection (e.g. used for centering UTM projection on faults",
    )
    parser.add_argument(
        "--dGrid",
        nargs=1,
        metavar=("dGrid"),
        default=[1.0],
        help="distance between consecutive parallel or meridians drawn",
        type=float,
    )
    parser.add_argument(
        "--beachSize",
        nargs=1,
        metavar=("dGrid"),
        default=[1.0],
        help="adjustement factor for beach ball sizes",
        type=float,
    )

    parser.add_argument(
        "--scalebarSize",
        nargs=1,
        metavar=("dGrid"),
        default=[100.0],
        help="size of scale bar in km",
        type=float,
    )
    parser.add_argument(
        "--shift",
        default=0.0,
        help="shift the rows of beach ball by shift for each row (useful if dip=90)",
        type=float,
    )
    parser.add_argument(
        "--MapBoundaries",
        nargs=4,
        metavar=("lonmin", "lonmax", "latmin", "latmax"),
        help="coordinates of map frame",
        type=float,
    )
    parser.add_argument(
        "--fault_edge",
        nargs=2,
        metavar=("fault_file", "projname"),
        help="add edge of fault to plot",
    )
    parser.add_argument(
        "--extension", nargs=1, default=(["svg"]), help="extension output file"
    )

    parser.add_argument(
        "--unicolor",
        dest="unicolor",
        action="store_true",
        help="beach balls all in blue",
    )

    # defer importing the heavy main function until execution
    def run(args):
        from seismic_waveform_factory.figure.plot_map_cmt import main

        main(args)

    parser.set_defaults(func=run)
