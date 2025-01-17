#!/usr/bin/env python3
from obspy import UTCDateTime, read
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
import pickle
import glob

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

source_files = config.get("GENERAL", "source_files").split(",")
source_files = [val.strip() for val in source_files]

seissol_outputs = config.get("GENERAL", "seissol_outputs", fallback=None)
if seissol_outputs:
    seissol_outputs = [val.strip() for val in seissol_outputs.split(",")]
    assert projection != None

all_files = []
for source_file in source_files:
    if os.path.isfile(source_file):
        if source_file.endswith(".h5"):
            all_files.append(source_file)
    # check for all point source files in the folder
    elif os.path.isdir(source_file):
        pattern = os.path.join(source_file, "PointSource*.h5")
        files_in_folder = glob.glob(pattern)
        print(files_in_folder)
        all_files.extend(files_in_folder)

source_files = sorted(list(set(all_files)))
print(f"{len(source_files)} source files found")

client_name = config.get("GENERAL", "client")
onset = config.get("GENERAL", "onset")
t1 = UTCDateTime(onset)
kind_vd = config.get("GENERAL", "kind")

software = config.get("GENERAL", "software").split(",")
# prevent result ['']
software = [item.strip() for item in software if item.strip()]

n_software = len(software)
prefix = f"plots/{setup_name}_{kind_vd}"
if "instaseis" in software:
    from instaseis_routines import generate_synthetics_instaseis

    db_name = config.get("GENERAL", "db")
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


hypo_lon = config.getfloat("GENERAL", "hypo_lon")
hypo_lat = config.getfloat("GENERAL", "hypo_lat")
scaling = config.getfloat("GENERAL", "scaling", fallback=1.0)
normalize = config.getboolean("GENERAL", "normalize", fallback=False)

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


n = len(source_files)
colors = extend_if_necessary(colors, n, "colors")
line_widths = extend_if_necessary(line_widths, n, "line_widths")


path_observations = config.get("GENERAL", "path_observations")
kind_misfit = config.get("GENERAL", "Misfit", fallback="rRMS")

plt.rcParams.update({"font.size": font_size})

os.makedirs(path_observations, exist_ok=True)

Pwave_tmin = -config.getfloat("P_WAVE", "t_before")
Pwave_tmax = config.getfloat("P_WAVE", "t_after")
Pwave_filter_fmin = 1.0 / config.getfloat("P_WAVE", "filter_tmax")
Pwave_filter_fmax = 1.0 / config.getfloat("P_WAVE", "filter_tmin")
Pwave_enabled = config.getboolean("P_WAVE", "enabled")
Pwave_ncol_per_component = config.getint("P_WAVE", "ncol_per_component")


SHwave_tmin = -config.getfloat("SH_WAVE", "t_before")
SHwave_tmax = config.getfloat("SH_WAVE", "t_after")
SHwave_filter_fmin = 1.0 / config.getfloat("SH_WAVE", "filter_tmax")
SHwave_filter_fmax = 1.0 / config.getfloat("SH_WAVE", "filter_tmin")
SHwave_enabled = config.getboolean("SH_WAVE", "enabled")
SHwave_ncol_per_component = config.getint("SH_WAVE", "ncol_per_component")


surface_waves_filter_fmin = 1.0 / config.getfloat("SURFACE_WAVES", "filter_tmax")
surface_waves_filter_fmax = 1.0 / config.getfloat("SURFACE_WAVES", "filter_tmin")
surface_waves_tmax = config.getfloat("SURFACE_WAVES", "tmax", fallback=None)
surface_waves_enabled = config.getboolean("SURFACE_WAVES", "enabled")
surface_waves_ncol_per_component = config.getint("SURFACE_WAVES", "ncol_per_component")
surface_waves_components = config.get("SURFACE_WAVES", "components").split(",")

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
n_kinematic_models = len(source_files)

Pwave = WaveformFigureGenerator(
    "P",
    Pwave_tmin,
    Pwave_tmax,
    Pwave_filter_fmin,
    Pwave_filter_fmax,
    Pwave_enabled,
    Pwave_ncol_per_component,
    nstations,
    ["Z"],
    n_software * n_kinematic_models,
    kind_misfit,
    colors,
    line_widths,
    scaling,
    normalize,
    relative_offset,
    annotations,
)
SHwave = WaveformFigureGenerator(
    "SH",
    SHwave_tmin,
    SHwave_tmax,
    SHwave_filter_fmin,
    SHwave_filter_fmax,
    SHwave_enabled,
    SHwave_ncol_per_component,
    nstations,
    ["T"],
    n_software * n_kinematic_models,
    kind_misfit,
    colors,
    line_widths,
    scaling,
    normalize,
    relative_offset,
    annotations,
)
surface_waves = WaveformFigureGenerator(
    "surface_waves",
    0.0,
    surface_waves_tmax,
    surface_waves_filter_fmin,
    surface_waves_filter_fmax,
    surface_waves_enabled,
    surface_waves_ncol_per_component,
    nstations,
    surface_waves_components,
    n_software * n_kinematic_models,
    kind_misfit,
    colors,
    line_widths,
    scaling,
    normalize,
    relative_offset,
    annotations,
)


components = ["E", "N", "Z"]

if Pwave_enabled and not (SHwave_enabled or surface_waves_enabled):
    # we do not need to compute E and N if Pwave only
    components = ["Z"]

gofall = [0 for i in range(len(source_files))]

list_synthetics_all = []
t_obs_before, t_obs_after = 100, 400

if seissol_outputs:
    t_obs_before, t_obs_after = 100, 400
    from seissol_receiver_processing import collect_seissol_synthetics

    list_synthetics = collect_seissol_synthetics(
        seissol_outputs, station_coords, projection, t1
    )
    list_synthetics_all += list_synthetics


if "axitra" in software and source_files:
    sys.path.append(path2axitra)
    from axitra_routines import generate_synthetics_axitra

    list_synthetics = generate_synthetics_axitra(
        source_files,
        station_coords,
        t1,
        kind_vd,
        fmax,
        duration,
        velocity_model_fname,
        path2axitra,
    )
    list_synthetics_all += list_synthetics
    t_obs_before, t_obs_after = 100, 400
    if not surface_waves_tmax:
        surface_waves.tmax = duration

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
    if not surface_waves_tmax:
        surface_waves.tmax = duration

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
    )
    list_synthetics_all += list_synthetics

    duration_synthetics = (
        list_synthetics[0][0].stats.endtime - list_synthetics[0][0].stats.starttime
    )

    t_obs_before, t_obs_after = 1000, duration_synthetics + 1000
    if not surface_waves_tmax:
        surface_waves.tmax = duration_synthetics

starttime = t1 - t_obs_before
endtime = t1 + t_obs_after
retrieved_waveforms = retrieve_waveforms(
    station_codes,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    processed_data,
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
    for st in [*lst, st_obs0]:
        for tr in st:
            tr.stats.back_azimuth = azimuth
            tr.stats.distance = dist

    if Pwave.enabled:
        tP = estimate_travel_time(hypo_depth_in_km, dist, station, "P")
        Pwave.add_plot_station(st_obs0, lst, t1 + tP, ins)

    if SHwave.enabled:
        tS = estimate_travel_time(hypo_depth_in_km, dist, station, "S")
        SHwave.add_plot_station(st_obs0, lst, t1 + tS, ins)

    if surface_waves.enabled:
        surface_waves.add_plot_station(st_obs0, lst, t1, ins)

print("goodness of fit (gof) per station:")
df_merged = merge_gof_dfs(Pwave, SHwave, surface_waves)

# Sort the column names alphabetically starting from "surface_waves_E0"
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
nsources = len(source_files)
df_station_average["source_file"] = [source_files[k % nsources] for k in file_id]
print(df_station_average)

fname = "gof_average.pkl"
df_station_average.to_pickle(fname)
print(f"done writing {fname}")

if not os.path.exists("plots"):
    os.makedirs("plots")

if Pwave.enabled:
    fname_P = f"{prefix}_Pwave.{ext}"
    Pwave.finalize_and_save_fig(fname_P)

if SHwave.enabled:
    fname_SH = f"{prefix}_SHwave.{ext}"
    SHwave.finalize_and_save_fig(fname_SH)

if surface_waves.enabled:
    fname_surface_waves = f"{prefix}_surface_waves_signal.{ext}"
    surface_waves.finalize_and_save_fig(fname_surface_waves)
