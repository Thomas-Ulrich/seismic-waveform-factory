#!/usr/bin/env python3
import glob
import argparse
import configparser
import os
import shutil


from waveform_figure_utils import compile_station_coords_main, estimate_travel_time
from obspy import UTCDateTime
from obspy.geodetics import locations2degrees
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

    station_codes = config.get("GENERAL", "stations").split(",")
    path_observations = config.get("GENERAL", "path_observations")

    print("updating strong_motion_waves.json with new durations")

    client_name = config.get("GENERAL", "client")
    onset = config.get("GENERAL", "onset")
    station_file = config.get("GENERAL", "station_file", fallback=None)
    hypo_lon = config.getfloat("GENERAL", "hypo_lon")
    hypo_lat = config.getfloat("GENERAL", "hypo_lat")
    hypo_depth_in_km = config.getfloat("GENERAL", "hypo_depth_in_km")

    t1 = UTCDateTime(onset)
    station_coords = compile_station_coords_main(
        station_codes, station_file, client_name, t1
    )
    station_2_full_code = {
        full_code.split(".")[1]: full_code for full_code in station_codes
    }
    print(station_2_full_code)
    with open(fn, "r") as f:
        data = json.load(f)
    # Update the duration for each entry
    for entry in data:
        station_name = entry["name"]
        station_code = station_2_full_code[station_name]
        lon, lat = station_coords[station_code]
        dist = locations2degrees(lat1=lat, long1=lon, lat2=hypo_lat, long2=hypo_lon)
        tP = estimate_travel_time(hypo_depth_in_km, dist, station, "P")
        if tP == 0.0:
            # station too close
            tP = 5.0
        entry["duration"] = int((tP + 45.0) / 2.0) * 10

    with open(fn, "w") as f:
        json.dump(data, f, indent=4)
    print(f"done updating {fn}")

