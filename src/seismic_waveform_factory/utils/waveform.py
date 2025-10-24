import functools as ft
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel

from seismic_waveform_factory.waveform.retrieve import get_station_data


def get_station_name_from_mseed(file_path):
    try:
        stream = read(file_path)
        # Assuming all traces in the stream belong to the same station
        station_name = stream[0].stats.station
        network_name = stream[0].stats.network
        return f"{network_name}.{station_name}"
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None


def get_station_files_dict(directory):
    station_files = {}
    for file_name in os.listdir(directory):
        if file_name.endswith(".mseed"):
            file_path = os.path.join(directory, file_name)
            station_name = get_station_name_from_mseed(file_path)
            if station_name:
                station_files[station_name] = file_path
    nfiles = len(station_files)
    print(f"found {nfiles} waveform files in {directory}")
    return station_files


def compile_station_coords_csv(station_codes, station_file):
    if os.path.exists(station_file):
        df = pd.read_csv(station_file)
        # Assert that all expected columns are present in the DataFrame
        expected_columns = ["station", "lon", "lat", "network"]
        for column in expected_columns:
            assert (
                column in df.columns
            ), f"Column '{column}' is missing from the DataFrame."

        return {
            f"{row['network']}.{row['station']}": (row["lon"], row["lat"])
            for _, row in df.iterrows()
        }


def extract_station_coords_from_dict(station_codes, station_coords_all, station_file):
    station_coords = {}
    for station in station_codes:
        if station in station_coords_all:
            station_coords[station] = station_coords_all[station]
        else:
            print(f"{station} not found in {station_file}")
    return station_coords


def compile_missing_stations(station_codes, station_coords_all):
    missing_stations = []
    for station in station_codes:
        if station not in station_coords_all:
            missing_stations += [station]
    return missing_stations


def download_station_coords(station_codes, client_name, t1):
    list_inventory = compile_list_inventories(client_name, station_codes, t1)
    print([f"{inv[0][0].code}" for inv in list_inventory])
    return compile_station_coords(list_inventory)


def compile_station_coords_main(station_codes, station_file, client_name, t1):
    station_coords = {}
    if station_file:
        station_coords_all = compile_station_coords_csv(station_codes, station_file)
        station_coords = extract_station_coords_from_dict(
            station_codes, station_coords_all, station_file
        )
    if len(station_coords) < len(station_codes):
        missing_stations = compile_missing_stations(station_codes, station_coords)
        print(missing_stations)
        downloaded_coords = download_station_coords(missing_stations, client_name, t1)
        station_coords = {**station_coords, **downloaded_coords}
    return station_coords


def parse_network_station(netStaCode):
    parts = netStaCode.split(".")
    return ("*", parts[0]) if len(parts) == 1 else (parts[0], parts[1])


def compile_list_inventories(client_name, station_codes, t1):
    # Prepare cache directory
    cache_dir = "observations"
    os.makedirs(cache_dir, exist_ok=True)

    def fetch_inventory(netStaCode):
        network, station = parse_network_station(netStaCode)
        # response because we will need it anyway at some point
        return get_station_data(
            client_name, network, [station], "*", "response", t1, t1 + 100.0, cache_dir
        )

    # Use ThreadPoolExecutor to process station codes in parallel
    list_inventory = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        # Submit tasks to the executor
        futures = {
            executor.submit(fetch_inventory, code): code for code in station_codes
        }

        # Collect results as they complete
        for future in as_completed(futures):
            netStaCode = futures[future]
            try:
                inventory = future.result()
                if inventory:
                    list_inventory.append(inventory)
            except Exception as e:
                print(f"Error processing {netStaCode}: {e}")

    return list_inventory


def compile_station_coords(list_inventory):
    station_coords = {}
    for ins, inv in enumerate(list_inventory):
        sta = inv[0][0]
        code = f"{inv[0].code}.{sta.code}"
        station_coords[code] = (sta.longitude, sta.latitude)
    return station_coords


def reorder_station_coords_from_azimuth(station_coords, hypo_lon, hypo_lat):
    # Reorder station based on azimuth
    d_azimuth = {}
    for code, lonlat in station_coords.items():
        longitude, latitude = lonlat
        d_azimuth[code] = gps2dist_azimuth(
            lat1=latitude, lon1=longitude, lat2=hypo_lat, lon2=hypo_lon
        )[2]
    ordered_azimuth = {
        k: v for k, v in sorted(d_azimuth.items(), key=lambda item: item[1])
    }
    reordered_station_coords = {}
    for key in ordered_azimuth:
        reordered_station_coords[key] = station_coords[key]
    return reordered_station_coords


def estimate_travel_time(source_depth_in_km, distance_in_degree, station, phase="P"):
    taupModel = "ak135"
    model = TauPyModel(model=taupModel)
    tP = model.get_travel_times(
        source_depth_in_km=source_depth_in_km,
        distance_in_degree=distance_in_degree,
        phase_list=[phase],
    )
    if not tP:
        print(f"no P wave at station {station}")
        tP = 0.0
    else:
        tP = tP[0].time
    return tP


def merge_gof_dfs(wf_plots):
    gofall_dfs = []
    for wf_plot in wf_plots:
        if wf_plot.enabled:
            df = wf_plot.gof_df.drop(columns=["distance", "azimuth"], errors="ignore")
            gofall_dfs.append(df)
    df_final = ft.reduce(
        lambda left, right: pd.merge(left, right, on="station", how="outer"), gofall_dfs
    )
    return df_final
