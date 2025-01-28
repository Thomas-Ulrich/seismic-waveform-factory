#!/usr/bin/env python3
import glob
import argparse
import configparser
import os
import shutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="move selected files to folder selected sac"
    )
    parser.add_argument("config_file", help="config file describing event and stations")
    args = parser.parse_args()

    config = configparser.ConfigParser()
    assert os.path.isfile(args.config_file), f"{args.config_file} not found"
    config.read(args.config_file)

    station_codes = config.get("GENERAL", "stations").split(",")
    path_observations = config.get("GENERAL", "path_observations")
    dest_folder = "selected_sac"
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
