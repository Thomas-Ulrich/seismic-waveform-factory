#!/usr/bin/env python3
import glob
import os
import sys
from itertools import cycle

import h5py
import matplotlib.pyplot as plt
import pandas as pd
from obspy import UTCDateTime
from obspy.geodetics import degrees2kilometers, locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
from seismic_waveform_factory.config.loader import ConfigLoader
from seismic_waveform_factory.config.schema import CONFIG_SCHEMA
from seismic_waveform_factory.config.utils import determine_config_scale
from seismic_waveform_factory.figure.generator import WaveformFigureGenerator
from seismic_waveform_factory.utils.waveform import (
    compile_station_coords_main,
    estimate_travel_time,
    get_station_files_dict,
    merge_gof_dfs,
    reorder_station_coords_from_azimuth,
)
from seismic_waveform_factory.waveform.retrieve import retrieve_waveforms


def main(args):
    # Ensure all rows and columns are displayed
    pd.set_option("display.max_rows", None)  # Show all rows
    pd.set_option("display.max_columns", None)  # Show all columns
    pd.set_option("display.width", 1000)  # Set a large width for the table
    pd.set_option(
        "display.colheader_justify", "center"
    )  # Center column headers (optional)

    assert os.path.isfile(args.config_file), f"{args.config_file} not found"
    cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)
    print(cfg["general"]["setup_name"])
    # print(cfg["waveform_plots"])

    def has_point_sources(filepath):
        with h5py.File(filepath, "r") as f:
            xyz_has_data = "xyz" in f and f["xyz"].shape[0] > 0
            return xyz_has_data

    def collect_synthetic_source_files(wf_syn_config):
        source_files = []
        all_files = []
        if "source_files" not in wf_syn_config.keys():
            return source_files
        for source_file in wf_syn_config["source_files"]:
            if os.path.isfile(source_file):
                if source_file.endswith(".h5") or source_file.endswith(".param"):
                    all_files.append(source_file)
            # check for all point source files in the folder
            elif os.path.isdir(source_file):
                pattern = os.path.join(source_file, "*ointSource*.h5")
                files_in_folder = glob.glob(pattern)
                print(files_in_folder)
                all_files.extend(files_in_folder)

        # remove duplicates but keep order
        seen = set()
        source_files = []
        for f in all_files:
            if f not in seen:
                seen.add(f)
                if f.endswith(".h5") and not has_point_sources(f):
                    print(f"removing {f} since it has no point sources")
                    continue
                source_files.append(f)

        syn_name = wf_syn_config["name"]
        print(f"{len(source_files)} source file(s) found for syn {syn_name}")
        return source_files

    syn_types = set()
    for wf_syn_config in cfg["synthetics"]:
        syn_types.add(wf_syn_config["type"])

    n_syn_model = {}
    for wf_syn_config in cfg["synthetics"]:
        source_files = collect_synthetic_source_files(wf_syn_config)
        wf_syn_config["source_files"] = source_files

        syn_name = wf_syn_config["name"]
        n_syn_model[syn_name] = (
            len(wf_syn_config["source_files"])
            if wf_syn_config["type"] != "seissol"
            else len(wf_syn_config["outputs"])
        )

    def extend_if_necessary(colors, n, name):
        ncolors = len(colors)
        if ncolors < n:
            print(f"enhancing line {name} as not enough specified")
            cycol = cycle(colors)
            for i in range(ncolors, n):
                colors.append(next(cycol))
        return colors

    print("n_syn_model:", n_syn_model)
    n_syn_model_max = 0
    for key in n_syn_model.keys():
        n_syn_model_max += n_syn_model[key]

    cfg["general"]["line_colors"] = extend_if_necessary(
        cfg["general"]["line_colors"], n_syn_model_max, "colors"
    )
    cfg["general"]["line_widths"] = extend_if_necessary(
        cfg["general"]["line_widths"], n_syn_model_max, "line_widths"
    )

    plt.rcParams.update({"font.size": cfg["general"]["font_size"]})

    wf_plots = []
    for plt_id, wf_plot_config in enumerate(cfg["waveform_plots"]):
        if wf_plot_config["annotations"]["distance_unit"] == "auto":
            unit = "degree" if "instaseis" in syn_types else "km"
            wf_plot_config["annotations"]["distance_unit"] = unit
        n_syn_models = 0
        for syn_name in wf_plot_config["synthetics"]:
            n_syn_models += n_syn_model[syn_name]
        wf_plot = WaveformFigureGenerator(
            cfg["general"], wf_plot_config, n_syn_models, plt_id
        )
        wf_plots.append(wf_plot)

    one_plot_enabled = any(wf_plot.enabled for wf_plot in wf_plots)
    if not one_plot_enabled:
        print("all plot disabled, exiting")
        return

    os.makedirs(cfg["general"]["path_observations"], exist_ok=True)

    processed_data = {}
    if cfg["processed_waveforms"]["directory"]:
        processed_data = cfg["processed_waveforms"]
        processed_data["station_files"] = get_station_files_dict(
            cfg["processed_waveforms"]["directory"]
        )

    station_file = cfg["general"]["station_file"]
    client_name = cfg["general"]["client"]
    hypo = cfg["general"]["hypocenter"]
    t1 = UTCDateTime(hypo["onset"])
    projection = cfg["general"]["projection"]
    path_observations = cfg["general"]["path_observations"]

    station_codes = set()
    for wf_plot in wf_plots:
        if wf_plot.enabled:
            for code in wf_plot.plt_cfg["stations"]:
                station_codes.add(code)
    station_codes = list(station_codes)

    station_coords = compile_station_coords_main(
        station_codes, station_file, client_name, t1, path_observations
    )
    print(station_coords)
    sta_infos = {}

    for ins, station_code in enumerate(station_coords):
        lon, lat = station_coords[station_code]
        network, station = station_code.split(".")
        lon2, lat2 = hypo["lon"], hypo["lat"]
        dist = locations2degrees(lat1=lat, long1=lon, lat2=lat2, long2=lon2)
        dist_km = degrees2kilometers(dist)
        azimuth = gps2dist_azimuth(lat1=lat, lon1=lon, lat2=lat2, lon2=lon2)[2]
        sta_infos[station_code] = {}
        sta_infos[station_code]["coords"] = station_coords[station_code]
        sta_infos[station_code]["dist"] = dist
        sta_infos[station_code]["dist_km"] = dist_km
        sta_infos[station_code]["azimuth"] = azimuth

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

    syn_info = {}
    for i, wf_syn_config in enumerate(cfg["synthetics"]):
        syn_name = wf_syn_config["name"]
        assert (
            syn_name not in syn_info.keys()
        ), f"duplicated name ({syn_name}) in synthetics"
        syn_info[syn_name] = {}
        for wf_plot in wf_plots:
            if wf_plot.enabled and syn_name in wf_plot.plt_cfg["synthetics"]:
                kind = wf_plot.plt_cfg["kind"]
                if kind not in syn_info[syn_name].keys():
                    syn_info[syn_name][kind] = {"stations": set()}
                for code in wf_plot.plt_cfg["stations"]:
                    syn_info[syn_name][kind]["stations"].add(code)

    def generate_synthetics(wf_syn_config, station_coords, syn_type):
        source_files = wf_syn_config["source_files"]
        if syn_type == "instaseis":
            from seismic_waveform_factory.simulation.instaseis import (
                generate_synthetics_instaseis,
            )

            list_synthetics = generate_synthetics_instaseis(
                wf_syn_config["db"],
                source_files,
                station_coords,
                t1,
                kind_vd,
                components,
                wf_syn_config["path_computed_synthetics"],
                projection,
                [wf_syn_config["mode"]],
            )

        elif syn_type == "seissol":
            from seismic_waveform_factory.simulation.seissol import (
                collect_seissol_synthetics,
            )

            assert projection is not None
            list_synthetics = collect_seissol_synthetics(
                wf_syn_config["outputs"], station_coords, projection, t1, kind_vd
            )

        elif syn_type == "axitra":
            for my_path in wf_syn_config["path"]:
                if os.path.isfile(os.path.join(my_path, "axitra")):
                    sys.path.append(my_path.strip())
                    break
            from seismic_waveform_factory.simulation.axitra import (
                generate_synthetics_axitra,
            )

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
            from seismic_waveform_factory.simulation.pyprop8 import (
                generate_synthetics_pyprop8,
            )

            list_synthetics = generate_synthetics_pyprop8(
                source_files,
                station_coords,
                t1,
                kind_vd,
                wf_syn_config["fmax"],
                wf_syn_config["duration"],
                wf_syn_config["velocity_model"],
            )
        else:
            raise ValueError("unknown synthetics type {syn_type}")
        try:
            duration = (
                list_synthetics[0][0].stats.endtime
                - list_synthetics[0][0].stats.starttime
            )
        except IndexError:
            # todo change hardcoded
            duration = 3600.0

        return list_synthetics, duration

    for i, wf_syn_config in enumerate(cfg["synthetics"]):
        syn_name = wf_syn_config["name"]
        for kind_vd in syn_info[syn_name].keys():
            selected_station_coords = {
                code: station_coords[code]
                for code in syn_info[syn_name][kind_vd]["stations"]
            }
            selected_station_coords = reorder_station_coords_from_azimuth(
                selected_station_coords, hypo["lon"], hypo["lat"]
            )

            syn_type = wf_syn_config["type"]
            st, duration = generate_synthetics(
                wf_syn_config, selected_station_coords, syn_type
            )
            syn_info[syn_name][kind_vd]["stream"] = st
            syn_info[syn_name][kind_vd]["duration"] = duration

    print(syn_info)

    merged = {}
    for syn_name, syn_types in syn_info.items():
        for kind, info in syn_types.items():
            if kind not in merged:
                merged[kind] = {"stations": set(), "duration": 0.0}
            merged[kind]["stations"].update(info["stations"])
            merged[kind]["duration"] = max(merged[kind]["duration"], info["duration"])

    config_scale = determine_config_scale(cfg)
    is_regional = config_scale["regional"]

    retrieved_waveforms = {}
    for kind_vd in merged.keys():
        duration = merged[kind_vd]["duration"]
        extra_time = max(100.0, 0.1 * duration)
        starttime = t1 - extra_time
        endtime = t1 + duration + extra_time
        retrieved_waveforms[kind_vd] = retrieve_waveforms(
            merged[kind_vd]["stations"],
            client_name,
            kind_vd,
            path_observations,
            starttime,
            endtime,
            processed_data=processed_data,
            is_regional=is_regional,
        )

    phase_dic = {"p": "P", "sh": "S"}
    for wf_plot in wf_plots:
        if wf_plot.enabled:
            print(wf_plot.plt_cfg)
            kind_vd = wf_plot.plt_cfg["kind"]
            selected_station_coords = {
                code: station_coords[code] for code in wf_plot.plt_cfg["stations"]
            }

            selected_station_coords = reorder_station_coords_from_azimuth(
                selected_station_coords, hypo["lon"], hypo["lat"]
            )

            for ins, code in enumerate(selected_station_coords.keys()):
                lst = []
                network, station = code.split(".")
                st_obs0 = retrieved_waveforms[kind_vd][code]
                for syn_name in wf_plot.plt_cfg["synthetics"]:
                    for st in syn_info[syn_name][kind_vd]["stream"]:
                        lst += [st.select(station=station)]

                for st in [*lst, st_obs0]:
                    for tr in st:
                        tr.stats.back_azimuth = sta_infos[code]["azimuth"]
                        tr.stats.distance = sta_infos[code]["dist"]
                        tr.stats.distance_km = sta_infos[code]["dist_km"]

                dist = sta_infos[code]["dist"]
                if wf_plot.plot_type in ["p", "sh"]:
                    t_phase = estimate_travel_time(
                        hypo["depth_in_km"], dist, station, phase_dic[wf_plot.plot_type]
                    )
                    wf_plot.set_estimated_travel_time(t_phase)
                else:
                    t_phase = 0.0
                wf_plot.add_plot_station(st_obs0, lst, t1 + t_phase, ins)

    src_loop_up = {}
    for wf_plot in wf_plots:
        if wf_plot.enabled:
            print(wf_plot.gof_df)
            src = []
            syn_names = wf_plot.plt_cfg["synthetics"]
            for wf_syn_config in cfg["synthetics"]:
                name = wf_syn_config["name"]
                if wf_syn_config["name"] in syn_names:
                    pt_sources = wf_syn_config.get("source_files") or wf_syn_config.get(
                        "outputs", []
                    )
                    src.extend([(name, pt_source) for pt_source in pt_sources])
            src_loop_up[f"{wf_plot.plt_id}"] = src
    print(src_loop_up)

    print("goodness of fit (gof) per station:")
    df_merged = merge_gof_dfs(wf_plots)

    if list(df_merged.columns) == ["station"]:
        print("not computing gof_per_station as df_merged is empty")
    else:
        fname = "gof_per_station.pkl"
        df_merged.to_pickle(fname)
        print(f"done writing {fname}")

        df_merged.drop(columns=["station"], inplace=True)

        df_station_average = df_merged.mean(axis=0).to_frame(name="gofa").reset_index()
        df_station_average = df_station_average.rename(columns={"index": "gofa_name"})
        file_id = (
            df_station_average["gofa_name"]
            .str.extract(r"(\d+)$", expand=False)
            .astype(int)
        )
        plot_id = (
            df_station_average["gofa_name"]
            .str.extract(r"(\d+)", expand=False)
            .astype(int)
        )
        df_station_average["plot_id"] = plot_id
        point_srcs = [
            src_loop_up[f"{p_id}"][f_id] for (p_id, f_id) in zip(plot_id, file_id)
        ]
        df_station_average["src"] = point_srcs

        df_station_average["sim_id"] = (
            df_station_average["src"]
            .apply(lambda x: x[1] if isinstance(x, tuple) and len(x) > 1 else str(x))
            .str.extract(r"dyn_(\d+)", expand=False)
            .fillna(-1)
            .astype(int)
        )

        print("average across all stations gof:")
        print(df_station_average)

        setup_name = cfg["general"]["setup_name"]
        fname = f"gof_{setup_name}_waveforms_average.pkl"
        df_station_average.to_pickle(fname)
        print(f"done writing {fname}")

    if not os.path.exists("plots"):
        os.makedirs("plots")

    for wf_plot in wf_plots:
        if wf_plot.enabled:
            wf_plot.finalize_and_save_fig()
