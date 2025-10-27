import argparse

import argcomplete

from seismic_waveform_factory.cmt.compute_multi_cmt_subparser import (
    add_parser as cmc_add_parser,
)
from seismic_waveform_factory.cmt.get_focal_mechanism_multi_cmt_subparser import (
    add_parser as gfm_add_parser,
)
from seismic_waveform_factory.cmt.plot_stf_cmt_subparser import (
    add_parser as psc_add_parser,
)
from seismic_waveform_factory.figure.generate_figure_synthetics_subparser import (
    add_parser as pws_add_parser,
)
from seismic_waveform_factory.figure.generate_legend_box_subparser import (
    add_parser as glb_add_parser,
)
from seismic_waveform_factory.figure.generate_map_with_waveforms_subparser import (
    add_parser as gmw_add_parser,
)
from seismic_waveform_factory.figure.plot_map_cmt_subparser import (
    add_parser as pmc_add_parser,
)
from seismic_waveform_factory.figure.plot_station_map_subparser import (
    add_parser as psm_add_parser,
)
from seismic_waveform_factory.geo.select_stations_subparser import (
    add_parser as ss_add_parser,
)

# done
from seismic_waveform_factory.seissol_utils.generate_station_file_subparser import (
    add_parser as gss_add_parser,
)

# not done


def main():
    parser = argparse.ArgumentParser(
        prog="swf", description="Seismic Waveform Factory CLI"
    )
    subparsers = parser.add_subparsers(title="subcommands", dest="command")
    subparsers.required = True

    # Register subcommands
    for add_sub in [
        cmc_add_parser,
        gfm_add_parser,
        pmc_add_parser,
        psc_add_parser,
        ss_add_parser,
        pws_add_parser,
        gss_add_parser,
        gmw_add_parser,
        psm_add_parser,
        glb_add_parser,
    ]:
        add_sub(subparsers)

    # Enable autocomplete
    argcomplete.autocomplete(parser)

    # Parse and dispatch
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
