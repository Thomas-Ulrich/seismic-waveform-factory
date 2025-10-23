#!/usr/bin/env python3
import argparse
import glob
import os

from obspy import read_inventory
from pyproj import Transformer
from seismic_waveform_factory.config.loader import ConfigLoader
from seismic_waveform_factory.config.schema import CONFIG_SCHEMA
from seismic_waveform_factory.config.utils import categorize_stations_by_scale

parser = argparse.ArgumentParser(description=("generate seissol station file"))
parser.add_argument("config_file", help="configuration file")
args = parser.parse_args()

assert os.path.isfile(args.config_file), f"{args.config_file} not found"
cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)
print(cfg["general"]["setup_name"])


# Extract station information
station_list = categorize_stations_by_scale(cfg)["regional"]
print(station_list)

projection = cfg["general"]["projection"]
# Set paths
waveform_folder = cfg["general"]["path_observations"]

# Initialize station dictionary
stations = {}

# Find waveform and response files
for station in station_list:
    net, sta = station.split(".")
    # Find response file
    response_file = glob.glob(f"{waveform_folder}{net}_{sta}_response.xml")
    if not response_file:
        print(f"Warning: No response found for {station}")
        continue
    response_file = response_file[0]

    # Store in dictionary
    stations[station] = {"response": response_file}

transformer = Transformer.from_crs(
    "epsg:4326",
    projection,
    always_xy=True,
)
coords = {}
coords_proj = {}
# Plot stations with waveforms
for station, info in stations.items():
    net, sta = station.split(".")
    inv = read_inventory(info["response"])
    y, x = (
        inv.networks[0].stations[0].latitude,
        inv.networks[0].stations[0].longitude,
    )
    coords[station] = [x, y]
    x, y = transformer.transform(x, y)
    coords_proj[station] = [x, y]

print(coords)
print(coords_proj)
fn = "seissol_receivers.txt"
with open(fn, "w") as fid:
    for station, coord in coords_proj.items():
        x, y = coord
        fid.write(f"{x} {y} {0}\n")

print(f"done writing {fn}")
