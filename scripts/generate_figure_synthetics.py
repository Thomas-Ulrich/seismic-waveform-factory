#!/usr/bin/env python3
import argparse
import configparser
import glob
import os
import sys
from itertools import cycle

import matplotlib.pyplot as plt
import pandas as pd
from config_loader import ConfigLoader
from config_schema import CONFIG_SCHEMA
from obspy import UTCDateTime
from obspy.geodetics import degrees2kilometers, locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
from retrieve_waveforms import retrieve_waveforms
from waveform_figure_generator import WaveformFigureGenerator
from waveform_figure_utils import (
    compile_station_coords_main,
    estimate_travel_time,
    get_station_files_dict,
    merge_gof_dfs,
    reorder_station_coords_from_azimuth,
)

# Ensure all rows and columns are displayed
pd.set_option("display.max_rows", None)  # Show all rows
pd.set_option("display.max_columns", None)  # Show all columns
pd.set_option("display.width", 1000)  # Set a large width for the table
pd.set_option("display.colheader_justify", "center")  # Center column headers (optional)


parser = argparse.ArgumentParser(
    description=("generate comparison plots of seismic waveform synthetics")
)
parser.add_argument("config_file", help="configuration file")
args = parser.parse_args()


assert os.path.isfile(args.config_file), f"{args.config_file} not found"
cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)
print(cfg["general"]["setup_name"])
# print(cfg["waveform_plots"])


def collect_synthetic_source_files(syn_config):
    source_files = []
    all_files = []
    for source_file in wf_syn_config["source_files"]:
        if os.path.isfile(source_file):
            if source_file.endswith(".h5"):
                all_files.append(source_file)
        # check for all point source files in the folder
        elif os.path.isdir(source_file):
            pattern = os.path.join(source_file, "*ointSource*.h5")
            files_in_folder = glob.glob(pattern)
            print(files_in_folder)
            all_files.extend(files_in_folder)
    source_files = sorted(list(set(all_files)))
    print(f"{len(source_files)} source file(s) found for syn {i} ({syn_type})")
    return source_files


syn_types = set()
n_syn_model = 0
for i, wf_syn_config in enumerate(cfg["synthetics"]):
    syn_type = wf_syn_config["type"]
    syn_types.add(syn_type)

    source_files = collect_synthetic_source_files(wf_syn_config)
    wf_syn_config["source_files"] = source_files

    if syn_type == "seissol":
        assert len(wf_syn_config["source_files"]) == 0
    if syn_type != "seissol":
        assert len(wf_syn_config["outputs"]) == 0
    n_syn_model += len(wf_syn_config["source_files"]) + len(wf_syn_config["outputs"])


print(syn_types)

if cfg["general"]["annotations"]["distance_unit"] == "auto":
    unit = "degree" if "instaseis" in syn_types else "km"
    cfg["general"]["annotations"]["distance_unit"] = unit


def extend_if_necessary(colors, n, name):
    ncolors = len(colors)
    if ncolors < n:
        print(f"enhancing line {name} as not enough specified")
        cycol = cycle(colors)
        for i in range(ncolors, n):
            colors.append(next(cycol))
    return colors


print("n_syn_model:", n_syn_model)

cfg["general"]["line_colors"] = extend_if_necessary(
    cfg["general"]["line_colors"], n_syn_model, "colors"
)
cfg["general"]["line_widths"] = extend_if_necessary(
    cfg["general"]["line_widths"], n_syn_model, "line_widths"
)

plt.rcParams.update({"font.size": cfg["general"]["font_size"]})


wf_plots = []
for wf_plot_config in cfg["waveform_plots"]:
    wf_plot = WaveformFigureGenerator(cfg["general"], wf_plot_config, n_syn_model)
    wf_plots.append(wf_plot)

os.makedirs(cfg["general"]["path_observations"], exist_ok=True)

processed_data = {}
"""
# To be reintegrated later
if config.has_section("PROCESSED_WAVEFORMS"):
    processed_data["directory"] = config.get("PROCESSED_WAVEFORMS", "directory")
    processed_data["wf_kind"] = config.get("PROCESSED_WAVEFORMS", "wf_kind")
    processed_data["wf_factor"] = config.getfloat(
        "PROCESSED_WAVEFORMS", "wf_factor", fallback=1.0
    )
    processed_data["station_files"] = get_station_files_dict(
        processed_data["directory"]
    )
"""


station_file = cfg["general"]["station_file"]
station_codes = cfg["general"]["stations"]
client_name = cfg["general"]["client"]
hypo = cfg["general"]["hypocenter"]
t1 = UTCDateTime(hypo["onset"])
projection = cfg["general"]["projection"]
kind_vd = cfg["general"]["kind"]
path_observations = cfg["general"]["path_observations"]

station_coords = compile_station_coords_main(
    station_codes, station_file, client_name, t1
)
station_coords = reorder_station_coords_from_azimuth(
    station_coords, hypo["lon"], hypo["lat"]
)
print(station_coords)


def collect_components(wf_plots):
    req = {c for wf in wf_plots for c in wf.components}
    components = []

    # Preserve E, N, Z order
    for v in ["E", "N", "Z"]:
        if v in req:
            components.append(v)

    # Special rule: if o, f, or T are requested, add E and N
    if {"o", "f", "T"} & req:
        for v in ["E", "N"]:
            if v not in components:
                components.append(v)

    return components


components = collect_components(wf_plots)

gofall = [0 for i in range(n_syn_model)]

list_synthetics_all = []
t_obs_before, t_obs_after = 100, 400


list_synthetics_all = []
for i, wf_syn_config in enumerate(cfg["synthetics"]):
    syn_type = wf_syn_config["type"]
    if syn_type == "instaseis":
        from instaseis_routines import generate_synthetics_instaseis

        list_synthetics = generate_synthetics_instaseis(
            wf_syn_config["db"],
            source_files,
            station_coords,
            t1,
            kind_vd,
            components,
            path_observations,
            projection,
            [wf_syn_config["mode"]],
        )

    elif syn_type == "seissol":
        from seissol_receiver_processing import collect_seissol_synthetics

        assert projection is not None
        list_synthetics = collect_seissol_synthetics(
            wf_syn_config["outputs"], station_coords, projection, t1, kind_vd
        )

    elif syn_type == "axitra":
        for my_path in wf_syn_config["path"]:
            if os.path.isfile(os.path.join(my_path, "axitra")):
                sys.path.append(my_path.strip())
                break
        from axitra_routines import generate_synthetics_axitra

        list_synthetics = generate_synthetics_axitra(
            source_files,
            station_coords,
            t1,
            kind_vd,
            wf_syn_config["fmax"],
            wf_syn_config["duration"],
            wf_syn_config["velocity_model"],
            my_path,
        )

    elif syn_type == "pyprop8":
        from pyprop8_routines import generate_synthetics_pyprop8

        list_synthetics = generate_synthetics_pyprop8(
            source_files,
            station_coords,
            t1,
            kind_vd,
            fmax,
            duration,
            velocity_model_fname,
        )
    else:
        raise ValueError("unknown synthetics type {syn_type}")
    list_synthetics_all += list_synthetics


"""
    #maybe not needed: axitra pyprop?
    t_obs_before, t_obs_after = 100, 400
    if not generic_wave.t_after:
        generic_wave.t_after = duration
    # instaseis 
    duration_synthetics = (
        list_synthetics[0][0].stats.endtime - list_synthetics[0][0].stats.starttime
    )

    t_obs_before, t_obs_after = 1000, duration_synthetics + 1000
    if not generic_wave.t_after:
        generic_wave.t_after = duration_synthetics

"""

starttime = t1 - t_obs_before
endtime = t1 + t_obs_after
retrieved_waveforms = retrieve_waveforms(
    station_codes,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    processed_data=processed_data,
)


for ins, station_code in enumerate(station_coords):
    lon, lat = station_coords[station_code]
    network, station = station_code.split(".")
    st_obs0 = retrieved_waveforms[f"{network}.{station}"]
    st_obs0.merge()

    network = st_obs0[0].stats.network

    lst = []
    for st_syn in list_synthetics_all:
        lst += [st_syn.select(station=station)]
    dist = locations2degrees(lat1=lat, long1=lon, lat2=hypo["lat"], long2=hypo["lon"])
    azimuth = gps2dist_azimuth(lat1=lat, lon1=lon, lat2=hypo["lat"], lon2=hypo["lon"])[
        2
    ]
    if cfg["general"]["annotations"]["distance_unit"] == "km":
        dist = degrees2kilometers(dist)
    for st in [*lst, st_obs0]:
        for tr in st:
            tr.stats.back_azimuth = azimuth
            tr.stats.distance = dist
            tr.stats.distance_unit = cfg["general"]["annotations"]["distance_unit"]

    phase_dic = {"p": "P", "sh": "S"}
    for wf_plot in wf_plots:
        if wf_plot.enabled:
            if wf_plot.signal_kind in ["p", "sh"]:
                t_phase = estimate_travel_time(
                    hypo["depth_in_km"], dist, station, phase_dic[wf_plot.signal_kind]
                )
                wf_plot.set_estimated_travel_time(t_phase)
            else:
                t_phase = 0.0
            wf_plot.add_plot_station(st_obs0, lst, t1 + t_phase, ins)


"""
print("goodness of fit (gof) per station:")
df_merged = merge_gof_dfs(Pwave, SHwave, generic_wave)

# Sort the column names alphabetically starting from "generic_wave_E0"
sorted_columns = sorted(df_merged.columns[df_merged.columns.get_loc("azimuth") + 1 :])
# Define the desired column order
desired_columns = ["station", "distance", "azimuth"] + sorted_columns

# Reorder the columns
df_merged = df_merged.reindex(columns=desired_columns)
print(df_merged)

fname = "gof_per_station.pkl"
df_merged.to_pickle(fname)
print(f"done writing {fname}")

df_merged.drop(columns=["station", "distance", "azimuth"], inplace=True)
print("station average gof:")
df_station_average = df_merged.mean(axis=0).to_frame(name="gofa").reset_index()
df_station_average = df_station_average.rename(columns={"index": "gofa_name"})
file_id = (
    df_station_average["gofa_name"]
    .str.extract(r"[^0-9]+([0-9]+)", expand=False)
    .astype(int)
)

source_files_inc_seissol = [seissol_outputs[k] for k in range(n_seissol_model)] + [
    source_files[(k - n_seissol_model) % n_kinematic_model]
    for k in range(n_seissol_model, n_syn_model)
]

pd.set_option("display.max_colwidth", None)
df_station_average["source_file"] = [source_files_inc_seissol[k] for k in file_id]
print(df_station_average.sort_values(by="gofa", ascending=False))

waveform_type = "teleseismic" if "instaseis" in software else "regional"
fname = f"gof_{waveform_type}_waveforms_average.pkl"
df_station_average.to_pickle(fname)
print(f"done writing {fname}")
"""

if not os.path.exists("plots"):
    os.makedirs("plots")

for wf_plot in wf_plots:
    if wf_plot.enabled:
        wf_plot.finalize_and_save_fig()
