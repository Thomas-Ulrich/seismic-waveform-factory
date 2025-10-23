#!/usr/bin/env python3
import argparse
import glob
import os
import shutil
from config.loader import ConfigLoader
from config.schema import CONFIG_SCHEMA

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "move selected files to folder selected sac and update "
            "strong_motion_waves.json if available"
        )
    )
    parser.add_argument("config_file", help="config file describing event and stations")
    parser.add_argument(
        "--output_folder", default=["selected_sac"], help="destination folder"
    )
    args = parser.parse_args()

    cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)

    station_codes = set()
    for plt_id, wf_plot_config in enumerate(cfg["waveform_plots"]):
        if wf_plot_config["enabled"]:
            for code in wf_plot_config["stations"]:
                station_codes.add(code)
    station_codes = list(station_codes)

    path_observations = cfg["general"]["path_observations"]

    dest_folder = args.output_folder
    os.makedirs(dest_folder, exist_ok=True)
    all_files = []
    for code in station_codes:
        net, station = code.split(".")
        patterns = [f"{code}*.sac", f"SAC_PZs_{net}_{station}*"]
        for pat in patterns:
            pattern = os.path.join(path_observations, pat)
            files_in_folder = glob.glob(pattern)
            all_files.extend(files_in_folder)

    print(all_files)

    for file in all_files:
        shutil.copy(file, dest_folder)
    print(f"done copying files to {dest_folder}")
