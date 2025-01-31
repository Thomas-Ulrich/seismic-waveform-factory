#!/usr/bin/env python3
import argparse
import configparser
import os
import shutil
from waveform_figure_utils import estimate_travel_time
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="update strong_motion_waves.json if available"
    )
    parser.add_argument("config_file", help="config file describing event and stations")
    args = parser.parse_args()

    config = configparser.ConfigParser()
    assert os.path.isfile(args.config_file), f"{args.config_file} not found"
    config.read(args.config_file)

    print("updating strong_motion_waves.json with new durations")
    hypo_depth_in_km = config.getfloat("GENERAL", "hypo_depth_in_km")

    fn = "strong_motion_waves.json"

    with open(fn, "r") as f:
        data = json.load(f)
    # Update the duration for each entry
    for entry in data:
        station_name = entry["name"]
        dist = entry["distance"]
        tP = estimate_travel_time(hypo_depth_in_km, dist, station_name, "P")
        if tP == 0.0:
            # station too close
            tP = 5.0
        entry["duration"] = int((tP + 50.0) / 2.0) * 10

    with open(fn, "w") as f:
        json.dump(data, f, indent=4)
    print(f"done updating {fn}")
