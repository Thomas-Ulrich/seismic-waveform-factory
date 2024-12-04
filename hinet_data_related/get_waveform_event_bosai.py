#!/usr/bin/env python3
import os
from datetime import timedelta
from HinetPy import Client
from obspy import UTCDateTime
import argparse

parser = argparse.ArgumentParser(
    description="download event data from bosai.jp using hinetpy"
)

parser.add_argument(
    "--origin_time_utc",
    nargs=1,
    required=True,
    metavar=("origin_time_of_eq"),
    help="time of the earthquake, e.g. 2024-01-01 07:10:09",
)
parser.add_argument(
    "--logging_data",
    nargs=2,
    metavar=("user", "password"),
    help="login data to Bosai",
    type=str,
    required=True,
)

parser.add_argument(
    "--hypocenter",
    nargs=2,
    metavar=("lon", "lat"),
    help="hypocenter location (center of the search radius)",
    required=True,
    type=float,
)
parser.add_argument(
    "--radius_range",
    nargs=2,
    metavar=("min", "max"),
    help="min and max radius to select stations",
    type=float,
    required=True,
)

args = parser.parse_args()

sorigin_time = args.origin_time_utc[0]
processed_date = sorigin_time.replace(" ", "_").replace(":","_")
outdir =  f"BOSAI_{processed_date}"

origin_time = UTCDateTime(sorigin_time).datetime
# UTC to Japan
origin_time += timedelta(hours=9)

username, password = args.logging_data
client = Client(username, password)
lon, lat = args.hypocenter
r1,r2 = args.radius_range
client.select_stations("0101", latitude=lat, longitude=lon, minradius=r1, maxradius=r2)

# skip if outdir already exists to avoid overwrite
if not os.path.exists(outdir):
    data, ctable = client.get_continuous_waveform("0101", origin_time-timedelta(minutes=3), 8, outdir=outdir)
