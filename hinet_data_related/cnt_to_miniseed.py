#!/usr/bin/env python3
from HinetPy import win32
import glob
from obspy import read, Stream
from obspy.signal.invsim import corn_freq_2_paz
from obspy.io.sac.sacpz import attach_paz
import os
import pandas as pd
import argparse

data = {"station": [], "lon": [], "lat": [], "network": []}
df = pd.DataFrame(data)

parser = argparse.ArgumentParser(description="convert kik net data to mseed")
parser.add_argument("data_folder", help="folder where the cnt and ch files are")
args = parser.parse_args()


folder_path = args.data_folder
cnt_files = [
    f"{folder_path}/{f}" for f in os.listdir(folder_path) if f.endswith(".cnt")
]
ch_files = [f"{folder_path}/{f}" for f in os.listdir(folder_path) if f.endswith(".ch")]
assert (len(cnt_files) == 1) and (len(ch_files) == 1)

sac_dir = f"{folder_path}_SAC"
win32.extract_sac(cnt_files[0], ch_files[0], outdir=sac_dir, with_sacpz=True)

# Collect all SAC files in the directory
sac_files = glob.glob(f"{sac_dir}/*.SAC")

# Initialize an empty stream
st = Stream()
paz_25hz = corn_freq_2_paz(25.0, damp=0.707)

# Read each SAC file and append to the stream
for sac_file in sac_files:
    st1 = read(sac_file)
    attach_paz(st1[0], f"{sac_file}_PZ")
    st1.simulate(paz_remove="self", paz_simulate=paz_25hz)
    # st1.differentiate()
    st += st1

stations = {f"{trace.stats.station}" for trace in st}
print(stations)
path_observations = f"{folder_path}_miniseed"

if not os.path.exists(path_observations):
    os.mkdir(path_observations)
t1 = st[0].stats.starttime

for sta in stations:
    selected_stream = st.select(station=sta)
    network = "BO"
    kind, station = sta.split(".")
    if kind == "N":
        # kind_vd = "acceleration"
        kind_vd = "velocity"
    else:
        kind_vd = "not_implemented"
    for tr in selected_stream:
        tr.stats.station = station
        tr.stats.network = network
        if tr.stats.channel.endswith("U"):
            tr.stats.channel = tr.stats.channel[:-1] + "Z"
        tr.data *= 1e-4
        tr.stats.starttime -= 9 * 3600  # 9 hours in seconds
    lat = tr.stats.sac.stla
    lon = tr.stats.sac.stlo
    # update for station file
    df.loc[len(df)] = {"station": station, "lon": lon, "lat": lat, "network": network}

    fname = f"{network}.{station}_{kind_vd}_{t1.date}.mseed"
    fullfname = os.path.join(path_observations, fname)

    selected_stream.write(fullfname, format="MSEED")

print(df)
fullfname = os.path.join(path_observations, "stations.csv")

df.to_csv(fullfname, index=False)
print(f"done writing {fullfname}")
