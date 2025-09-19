#!/usr/bin/env python3
from obspy import UTCDateTime
from obspy.geodetics import locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
import argparse
import configparser
import matplotlib.pyplot as plt
import os
from retrieve_waveforms import retrieve_waveforms
from waveform_figure_generator import WaveformFigureGenerator
import sys
from waveform_figure_utils import (
    compile_station_coords_main,
    estimate_travel_time,
    get_station_files_dict,
    merge_gof_dfs,
    reorder_station_coords_from_azimuth,
)
from itertools import cycle
import glob
import pandas as pd
from obspy.geodetics import degrees2kilometers

# Ensure all rows and columns are displayed
pd.set_option("display.max_rows", None)  # Show all rows
pd.set_option("display.max_columns", None)  # Show all columns
pd.set_option("display.width", 1000)  # Set a large width for the table
pd.set_option("display.colheader_justify", "center")  # Center column headers (optional)


parser = argparse.ArgumentParser(
    description=(
        "generate synthetics with instaseis or axitra given sources and stations"
    )
)
parser.add_argument("config_file", help="config file describing event and stations")
args = parser.parse_args()


config = configparser.ConfigParser()
assert os.path.isfile(args.config_file), f"{args.config_file} not found"
config.read(args.config_file)

setup_name = config.get("GENERAL", "setup_name")
ext = config.get("GENERAL", "figure_extension")
font_size = config.get("GENERAL", "font_size")
relative_offset = config.getfloat("GENERAL", "relative_offset", fallback=0.0)
annotations = config.get(
    "GENERAL", "annotations", fallback="distance,azimuth,misfit"
).split(",")

projection = config.get("GENERAL", "projection", fallback=None)
distance_unit = config.get("GENERAL", "distance_unit", fallback=None)

source_files = config.get("GENERAL", "source_files").split(",")
source_files = [val.strip() for val in source_files]

seissol_outputs = config.get("GENERAL", "seissol_outputs", fallback=[])
if seissol_outputs:
    seissol_outputs = [val.strip() for val in seissol_outputs.split(",")]
    assert projection is not None

all_files = []
for source_file in source_files:
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
print(f"{len(source_files)} source files found")

client_name = config.get("GENERAL", "client")
onset = config.get("GENERAL", "onset")
t1 = UTCDateTime(onset)
kind_vd = config.get("GENERAL", "kind")
global_legend_labels = config.get("GENERAL", "global_legend_labels", fallback=None)
if global_legend_labels:
    global_legend_labels = [item.strip() for item in global_legend_labels.split(";")]


software = config.get("GENERAL", "software").split(",")
# prevent result ['']
software = [item.strip() for item in software if item.strip()]

n_software = len(software)
prefix = f"plots/{setup_name}_{kind_vd}"
if "instaseis" in software:
    from instaseis_routines import generate_synthetics_instaseis

    db_name = config.get("GENERAL", "db")
    instaseis_modes = config.get("GENERAL", "instaseis_mode", fallback="classical")
    instaseis_modes = [item.strip() for item in instaseis_modes.split(",")]

    prefix = f"{prefix}_instaseis"

if ("axitra" in software) or ("pyprop8" in software):
    fmax = config.getfloat("GENERAL", "axitra_pyprop8_fmax")
    duration = config.getfloat("GENERAL", "axitra_pyprop8_duration")
    velocity_model_fname = config.get("GENERAL", "axitra_pyprop8_velocity_model")
    if "axitra" in software:
        path2axitra = config.get("GENERAL", "axitra_path")
    prefix = f"{prefix}_{fmax}Hz_{duration}s"

if seissol_outputs:
    prefix = f"{prefix}_{kind_vd}_seissol"

if not distance_unit:
    if ("axitra" in software) or ("pyprop8" in software) or seissol_outputs:
        distance_unit = "km"
    else:
        distance_unit = "degree"


hypo_lon = config.getfloat("GENERAL", "hypo_lon")
hypo_lat = config.getfloat("GENERAL", "hypo_lat")
scaling = config.getfloat("GENERAL", "scaling", fallback=1.0)
normalize = config.getboolean("GENERAL", "normalize", fallback=False)
shift_match_correlation = config.getboolean(
    "GENERAL", "shift_match_correlation", fallback=False
)

hypo_depth_in_km = config.getfloat("GENERAL", "hypo_depth_in_km")
station_codes_list = config.get("GENERAL", "stations")
station_codes = [v.strip() for v in station_codes_list.split(",")]
colors = config.get("GENERAL", "line_colors").split(",")
line_widths = [float(v) for v in config.get("GENERAL", "line_widths").split(",")]


def extend_if_necessary(colors, n, name):
    ncolors = len(colors)
    if ncolors < n:
        print(f"enhancing line {name} as not enough specified")
        cycol = cycle(colors)
        for i in range(ncolors, n):
            colors.append(next(cycol))
    return colors


n_seissol_model = len(seissol_outputs)
n_kinematic_model = len(source_files)
n_syn_model = n_software * n_kinematic_model + n_seissol_model
print("n_syn_model:", n_syn_model)

colors = extend_if_necessary(colors, n_syn_model, "colors")
line_widths = extend_if_necessary(line_widths, n_syn_model, "line_widths")


path_observations = config.get("GENERAL", "path_observations")
kind_misfit = config.get("GENERAL", "Misfit", fallback="min_shifted_normalized_rms")

valid_misfits = {
    "min_shifted_normalized_rms",
    "normalized_rms",
    "cross-correlation",
    "time-frequency",
}
if kind_misfit not in valid_misfits:
    raise ValueError(
        f"Invalid misfit kind: {kind_misfit}. Must be one of {valid_misfits}"
    )


plt.rcParams.update({"font.size": font_size})

os.makedirs(path_observations, exist_ok=True)


def get_wave_config(config, section_names):
    # determine section name from all possibilities
    section = None
    for section_name in section_names:
        if config.has_section(section_name):
            section = section_name
            break
    if not section:
        raise ValueError(f"none of {section_names} found in config")
    t_before = 0
    t_after = None
    if config.has_option(section, "t_before"):
        t_before = -config.getfloat(section, "t_before")
        t_after = config.getfloat(section, "t_after")
    elif config.has_option(section, "tmax"):
        t_before = config.getfloat(section, "tmin", fallback=0.0)
        t_after = config.getfloat(section, "tmax", fallback=None)

    taper = config.getboolean(section, "taper", fallback=True)
    filter_fmin = 1.0 / config.getfloat(section, "filter_tmax")
    filter_fmax = 1.0 / config.getfloat(section, "filter_tmin")
    fault_strike = config.getfloat(section, "fault_strike", fallback=None)
    enabled = config.getboolean(section, "enabled")
    ncol = config.getint(section, "ncol_per_component")

    # Optionally get components, fallback empty list if missing
    if config.has_option(section, "components"):
        components = [c.strip() for c in config.get(section, "components").split(",")]
        if "o" in components or "f" in components:
            assert fault_strike is not None, (
                "fault_strike parameter required for components='f' or 'o'"
                "(fault-parallel and normal)"
            )
    elif section in ["P_WAVE"]:
        components = ["Z"]
    elif section in ["SH_WAVE"]:
        components = ["T"]

    return {
        "t_before": t_before,
        "t_after": t_after,
        "taper": taper,
        "filter_fmin": filter_fmin,
        "filter_fmax": filter_fmax,
        "enabled": enabled,
        "ncol_per_component": ncol,
        "components": components,
        "fault_strike": fault_strike,
    }


p_window_config = get_wave_config(config, ["P_WAVE"])
s_window_config = get_wave_config(config, ["SH_WAVE"])
origin_time_window_config = get_wave_config(config, ["GENERIC_WAVE", "SURFACE_WAVES"])

processed_data = {}
if config.has_section("PROCESSED_WAVEFORMS"):
    processed_data["directory"] = config.get("PROCESSED_WAVEFORMS", "directory")
    processed_data["wf_kind"] = config.get("PROCESSED_WAVEFORMS", "wf_kind")
    processed_data["wf_factor"] = config.getfloat(
        "PROCESSED_WAVEFORMS", "wf_factor", fallback=1.0
    )
    processed_data["station_files"] = get_station_files_dict(
        processed_data["directory"]
    )


station_file = config.get("GENERAL", "station_file", fallback=None)
station_coords = compile_station_coords_main(
    station_codes, station_file, client_name, t1
)
station_coords = reorder_station_coords_from_azimuth(station_coords, hypo_lon, hypo_lat)
print(station_coords)

nstations = len(station_coords)


Pwave = WaveformFigureGenerator(
    "P",
    nstations,
    n_syn_model,
    kind_misfit,
    colors,
    line_widths,
    scaling,
    normalize,
    shift_match_correlation,
    relative_offset,
    annotations,
    global_legend_labels,
    **p_window_config,
)

SHwave = WaveformFigureGenerator(
    "SH",
    nstations,
    n_syn_model,
    kind_misfit,
    colors,
    line_widths,
    scaling,
    normalize,
    shift_match_correlation,
    relative_offset,
    annotations,
    global_legend_labels,
    **s_window_config,
)
generic_wave = WaveformFigureGenerator(
    "generic",
    nstations,
    n_syn_model,
    kind_misfit,
    colors,
    line_widths,
    scaling,
    normalize,
    shift_match_correlation,
    relative_offset,
    annotations,
    global_legend_labels,
    **origin_time_window_config,
)

components = ["E", "N", "Z"]

if Pwave.enabled and not (SHwave.enabled or generic_wave.enabled):
    # we do not need to compute E and N if Pwave only
    components = ["Z"]

gofall = [0 for i in range(n_syn_model)]

list_synthetics_all = []
t_obs_before, t_obs_after = 100, 400

if seissol_outputs:
    t_obs_before, t_obs_after = 100, 400
    from seissol_receiver_processing import collect_seissol_synthetics

    list_synthetics = collect_seissol_synthetics(
        seissol_outputs, station_coords, projection, t1, kind_vd
    )
    list_synthetics_all += list_synthetics


if "axitra" in software and source_files:
    for my_path in path2axitra.split(","):
        if os.path.isfile(os.path.join(my_path, "axitra")):
            sys.path.append(my_path.strip())
            break
    from axitra_routines import generate_synthetics_axitra

    list_synthetics = generate_synthetics_axitra(
        source_files,
        station_coords,
        t1,
        kind_vd,
        fmax,
        duration,
        velocity_model_fname,
        my_path,
    )
    list_synthetics_all += list_synthetics
    t_obs_before, t_obs_after = 100, 400
    if not generic_wave.t_after:
        generic_wave.t_after = duration

if "pyprop8" in software and source_files:
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
    list_synthetics_all += list_synthetics
    t_obs_before, t_obs_after = 100, 400
    if not generic_wave.t_after:
        generic_wave.t_after = duration

if "instaseis" in software and source_files:
    list_synthetics = generate_synthetics_instaseis(
        db_name,
        source_files,
        station_coords,
        t1,
        kind_vd,
        components,
        path_observations,
        projection,
        instaseis_modes,
    )
    list_synthetics_all += list_synthetics

    duration_synthetics = (
        list_synthetics[0][0].stats.endtime - list_synthetics[0][0].stats.starttime
    )

    t_obs_before, t_obs_after = 1000, duration_synthetics + 1000
    if not generic_wave.t_after:
        generic_wave.t_after = duration_synthetics

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
    dist = locations2degrees(lat1=lat, long1=lon, lat2=hypo_lat, long2=hypo_lon)
    azimuth = gps2dist_azimuth(lat1=lat, lon1=lon, lat2=hypo_lat, lon2=hypo_lon)[2]
    if distance_unit == "km":
        dist = degrees2kilometers(dist)
    for st in [*lst, st_obs0]:
        for tr in st:
            tr.stats.back_azimuth = azimuth
            tr.stats.distance = dist
            tr.stats.distance_unit = distance_unit

    if Pwave.enabled:
        tP = estimate_travel_time(hypo_depth_in_km, dist, station, "P")
        Pwave.set_estimated_travel_time(tP)
        Pwave.add_plot_station(st_obs0, lst, t1 + tP, ins)

    if SHwave.enabled:
        tS = estimate_travel_time(hypo_depth_in_km, dist, station, "S")
        SHwave.set_estimated_travel_time(tS)
        SHwave.add_plot_station(st_obs0, lst, t1 + tS, ins)

    if generic_wave.enabled:
        generic_wave.add_plot_station(st_obs0, lst, t1, ins)

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

if not os.path.exists("plots"):
    os.makedirs("plots")

if Pwave.enabled:
    fname_P = f"{prefix}_P.{ext}"
    Pwave.finalize_and_save_fig(fname_P)

if SHwave.enabled:
    fname_SH = f"{prefix}_SH.{ext}"
    SHwave.finalize_and_save_fig(fname_SH)

if generic_wave.enabled:
    fname_generic_wave = f"{prefix}_generic.{ext}"
    generic_wave.finalize_and_save_fig(fname_generic_wave)
