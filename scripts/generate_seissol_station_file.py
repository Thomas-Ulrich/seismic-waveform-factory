import numpy as np
import configparser
import os
import glob
import matplotlib.pyplot as plt
from obspy import read, read_inventory
from obspy.clients.nrl import NRL
from obspy.core import UTCDateTime
from pyproj import Transformer

# Load the configuration file
config = configparser.ConfigParser()
config.read("waveforms_config.ini")

# Extract station information
station_list = [st.strip() for st in config["GENERAL"]["stations"].split(",")]
projection = config["GENERAL"]["projection"]
# Set paths
waveform_folder = "observations/"

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
