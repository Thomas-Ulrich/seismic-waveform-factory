#!/usr/bin/env python3
from obspy import UTCDateTime
from obspy.clients.fdsn import Client, RoutingClient
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.inventory import Inventory
from obspy import read_inventory
import pandas as pd
import numpy as np
import argparse
import configparser
import os
import jinja2
import geopandas as gpd
from geopy.distance import geodesic
from shapely.geometry import Point
from shapely.ops import nearest_points
from pyproj import Transformer
from retrieve_waveforms import retrieve_waveforms_including_preprocessed
from plot_station_map import generate_station_map
from fault_processing import compute_shapely_polygon, get_fault_slip_coords
from scipy import spatial
from waveform_figure_utils import get_station_files_dict

np.random.seed(42)


def generate_station_df(inv0):
    df = pd.DataFrame(columns=["network", "station", "longitude", "latitude"])
    for i, net in enumerate(inv0):
        for j, sta in enumerate(net):
            new_row = {
                "network": net.code,
                "station": sta.code,
                "longitude": sta.longitude,
                "latitude": sta.latitude,
            }
            df.loc[len(df)] = new_row
    return df


def compute_dict_network_station(df):
    networks = set(df["network"])
    network_station = {}
    for net in networks:
        network_station[net] = []
    for index, row in df.iterrows():
        network_station[row["network"]].append(row["station"])
    return network_station


# Define a function to calculate the distance between two points using their latitude and longitude
def haversine_distance(lat1, lon1, lat2, lon2):
    coords_1 = (lat1, lon1)
    coords_2 = (lat2, lon2)
    return geodesic(coords_1, coords_2).km


def remove_synthetics_from_inventory(original_inv):
    new_inv = Inventory()

    for network in original_inv:
        # Create a new network
        new_network = network.copy()
        new_network.stations = []  # Clear the stations list
        for station in network:
            if not station.site.name:
                new_network.stations.append(station)
            elif "synthetic" not in station.site.name.lower():
                new_network.stations.append(station)

        # If the network has any stations left, add it to the new inventory
        if new_network.stations:
            new_inv.networks.append(new_network)
    return new_inv


def generate_geopanda_dataframe(df_stations, fault_info):
    df = pd.DataFrame(
        columns=["network", "station", "longitude", "latitude", "distance_km"]
    )
    if fault_info:
        projection = fault_info["projection"]
        transformer = Transformer.from_crs("epsg:4326", projection, always_xy=True)
        coords = fault_info["fault_slip_coords"] / 1e3
        tree = spatial.KDTree(coords)

    for _, row in df_stations.iterrows():
        if fault_info:
            x1, y1 = transformer.transform(row.longitude, row.latitude)
            min_distance, _ = tree.query([x1 / 1e3, y1 / 1e3, 0])
        else:
            min_distance = -1.0
        new_row = {
            "network": row.network,
            "station": row.station,
            "longitude": row.longitude,
            "latitude": row.latitude,
            "distance_km": min_distance,
        }
        df.loc[len(df)] = new_row
    df["code"] = df["network"] + "." + df["station"]
    df.sort_values(by="distance_km", inplace=True)
    df = df.reset_index(drop=True)
    # Convert the latitude and longitude columns to a GeoDataFrame
    geometry = gpd.points_from_xy(df.longitude, df.latitude)
    crs = "epsg:4326"
    return gpd.GeoDataFrame(df, crs=crs, geometry=geometry)


def select_closest_stations(available_stations, selected_stations, nstations):
    if len(selected_stations) == 0:
        print("selected_stations is empty, selcting the n closest nstation")
        selected_stations = available_stations.head(nstations)
        available_stations = available_stations.drop(selected_stations.index)

    nselected = len(selected_stations)
    remaining_stations = nstations - nselected

    return selected_stations.sort_index(), available_stations.sort_index()


def select_stations_most_distant(available_stations, selected_stations, nstations):
    if len(selected_stations) == 0:
        print("selected_stations is empty, selecting randomly from available_stations")
        selected_stations = available_stations.sample(1)
        available_stations = available_stations.drop(selected_stations.index)

    nselected = len(selected_stations)
    remaining_stations = nstations - nselected

    if remaining_stations <= 0:
        return selected_stations.sort_index(), available_stations.sort_index()

    for _ in range(remaining_stations):
        if available_stations.empty:
            print(
                f"GeoDataFrame is empty, {remaining_stations} stations will be missing"
            )
            break

        distances = []
        for _, row in available_stations.iterrows():
            distance = min(
                haversine_distance(
                    row.latitude, row.longitude, station.latitude, station.longitude
                )
                for station in selected_stations.itertuples(index=False)
            )
            distances.append((distance, row))

        furthest_station = max(distances, key=lambda x: x[0])[1]
        selected_stations = pd.concat(
            [selected_stations, furthest_station.to_frame().T]
        )
        available_stations = available_stations.drop(furthest_station.name)

    return selected_stations.sort_index(), available_stations.sort_index()


def select_telseismic_stations_aiming_for_azimuthal_coverage(inv0, event, nstations):
    df = pd.DataFrame(
        columns=[
            "network",
            "station",
            "longitude",
            "latitude",
            "distance",
            "backazimuth",
        ]
    )
    for i, net in enumerate(inv0):
        for j, sta in enumerate(net):
            distance, azimuth, backazimuth = gps2dist_azimuth(
                lat1=sta.latitude,
                lon1=sta.longitude,
                lat2=event["latitude"],
                lon2=event["longitude"],
            )
            new_row = {
                "network": net.code,
                "station": sta.code,
                "longitude": sta.longitude,
                "latitude": sta.latitude,
                "distance": distance,
                "backazimuth": backazimuth,
            }
            df.loc[len(df)] = new_row
    df.sort_values(by="backazimuth", inplace=True)
    df["ba_difference"] = df["backazimuth"].diff()
    df["code"] = df["network"] + "." + df["station"]
    df["ba_difference"] = df["ba_difference"].fillna(0)
    df["cumulative_ba"] = df["ba_difference"].cumsum()
    df = df.sort_values(by="cumulative_ba")
    df = df.reset_index(drop=True)
    # Values to pick
    target_values = np.linspace(0, 360, nstations + 1)[:-1]
    print(target_values)
    # Find te indices of the nearest values
    nearest_indices = []
    for val in target_values:
        nearest_indices.append((df["cumulative_ba"] - val).abs().idxmin())
    # Slice the DataFrame using the nearest indices
    df_selected = df.loc[nearest_indices]
    return df_selected, df


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Select stations ensuring optimal coverage for an earthquake."
    )
    parser.add_argument("config_file", help="Config file describing the event.")
    parser.add_argument(
        "number_stations", type=int, help="Number of stations to select (in total)"
    )
    parser.add_argument(
        "closest_stations", type=int, help="Number of the closest station to select."
    )

    parser.add_argument(
        "--azimuthal",
        action="store_true",
        help="Select stations based on back-azimuth (else based on distance).",
    )
    return parser.parse_args()


def read_config(config_file):
    config = configparser.ConfigParser()
    assert os.path.isfile(config_file), f"{config_file} not found"
    config.read(config_file)
    return config


def load_or_create_inventory(
    client, client_name, event, path_observations, spatial_range, t_before, t_after
):
    fn_inventory = f"{path_observations}/inv_{client_name}.xml"
    if os.path.exists(fn_inventory):
        inventory = read_inventory(fn_inventory)
    else:
        starttime = event["onset"] - t_before
        endtime = event["onset"] + t_after
        kargs = {}
        r0, r1 = spatial_range["radius"]
        if "source_lon_range" in spatial_range.keys():
            kargs["minlongitude"] = spatial_range["source_lon_range"][0] - r1
            kargs["maxlongitude"] = spatial_range["source_lon_range"][1] + r1
            kargs["minlatitude"] = spatial_range["source_lat_range"][0] - r1
            kargs["maxlatitude"] = spatial_range["source_lat_range"][1] + r1
        else:
            kargs["minradius"] = r0
            kargs["maxradius"] = r1
            kargs["latitude"] = event["latitude"]
            kargs["longitude"] = event["longitude"]
        print(kargs)
        inventory = client.get_stations(
            starttime=starttime,
            endtime=endtime,
            level="channel",
            channel="*Z",
            **kargs,
        )
        inventory.write(fn_inventory, format="STATIONXML")
    return inventory


def compute_min_max_coords(fault_info):
    # Initialize the transformer with the provided projection
    projection = fault_info["projection"]
    transformer = Transformer.from_crs(projection, "epsg:4326", always_xy=True)
    polygons = fault_info["polygons"]

    # Initialize min and max coordinates
    min_x, min_y = float("inf"), float("inf")
    max_x, max_y = float("-inf"), float("-inf")

    # Transform each polygon and compute min/max coordinates
    for poly in polygons:
        x1, y1 = poly.exterior.xy
        x1, y1 = transformer.transform(x1, y1)

        min_x = min(min_x, min(x1))
        min_y = min(min_y, min(y1))
        max_x = max(max_x, max(x1))
        max_y = max(max_y, max(y1))

    # Return the min and max coordinates
    return (min_x, max_x), (min_y, max_y)


if __name__ == "__main__":
    args = parse_arguments()
    config = read_config(args.config_file)
    setup_name = config.get("GENERAL", "setup_name")
    client_name = config.get("GENERAL", "client")
    onset = config.get("GENERAL", "onset")
    hypo_lon = config.getfloat("GENERAL", "hypo_lon")
    hypo_lat = config.getfloat("GENERAL", "hypo_lat")
    config_stations = config.get("GENERAL", "stations", fallback="")
    path_observations = config.get("GENERAL", "path_observations")
    kind_vd = config.get("GENERAL", "kind")
    software = config.get("GENERAL", "software").split(",")
    is_teleseismic = "axitra" not in software

    faultfname = config.get("GENERAL", "fault_filename", fallback=None)
    station_file = config.get("GENERAL", "station_file", fallback=None)

    processed_data = {}
    processed_data["directory"] = config.get(
        "GENERAL", "processed_waveforms", fallback=None
    )

    if processed_data["directory"]:
        processed_data["wf_kind"] = config.get("GENERAL", "processed_waveforms_kind")
        processed_data["wf_factor"] = config.getfloat(
            "GENERAL", "processed_waveforms_factor", fallback=1.0
        )
        processed_data["station_files"] = get_station_files_dict(
            processed_data["directory"]
        )

    spatial_range = {}

    if faultfname:
        polygons = compute_shapely_polygon(faultfname)
        fault_slip_coords = get_fault_slip_coords(faultfname)
        projection = config.get("GENERAL", "projection")
        fault_info = {}
        fault_info["projection"] = projection
        fault_info["fault_slip_coords"] = fault_slip_coords
        fault_info["polygons"] = polygons
        lon_range, lat_range = compute_min_max_coords(fault_info)
        spatial_range["source_lon_range"] = lon_range
        spatial_range["source_lat_range"] = lat_range

    else:
        fault_info = None

    if not is_teleseismic:
        duration = config.getfloat("GENERAL", "axitra_pyprop8_duration")
        t_before = 100
        t_after = duration + 100
        spatial_range["radius"] = [0.0, 2.5]
    else:
        import instaseis

        db_name = config.get("GENERAL", "db")
        db = instaseis.open_db(db_name)
        t_before = 1000
        t_after = db.info.length + 1000
        spatial_range["radius"] = [30, 90]

    if client_name in ["eida-routing", "iris-federator"]:
        client = RoutingClient(client_name)
    else:
        client = Client(client_name)

    # Define the event parameters
    event_time = UTCDateTime(onset)
    event = {"longitude": hypo_lon, "latitude": hypo_lat, "onset": event_time}
    starttime = event["onset"] - t_before
    endtime = event["onset"] + t_after

    os.makedirs(path_observations, exist_ok=True)

    if station_file:
        station_df = pd.read_csv(station_file)
        station_df.rename(columns={"lon": "longitude", "lat": "latitude"}, inplace=True)
        print(station_df)
    else:
        inventory = load_or_create_inventory(
            client,
            client_name,
            event,
            path_observations,
            spatial_range,
            t_before,
            t_before,
        )
        inventory = remove_synthetics_from_inventory(inventory)
        print("remove AM network")
        filtered_networks = [net for net in inventory.networks if net.code != "AM"]
        inventory = Inventory(networks=filtered_networks, source=inventory.source)
        station_df = generate_station_df(inventory)

    available_stations = generate_geopanda_dataframe(station_df, fault_info)
    projection = config.get("GENERAL", "projection", fallback="")
    if projection:
        # 60 closest stations are written for seissol output
        closest_stations = available_stations.head(60)
        transformer = Transformer.from_crs("epsg:4326", projection, always_xy=True)
        x1, y1 = transformer.transform(
            closest_stations["longitude"], closest_stations["latitude"]
        )
        n_seissol_station = len(closest_stations)
        for k in [n_seissol_station, 5]:
            fname = f"tmp/seissol_station_{k}.txt"
            os.makedirs("tmp", exist_ok=True)
            with open(fname, "w") as fid:
                for i in range(k):
                    fid.write(f"{x1[i]} {y1[i]} 0\n")
            print(f"done writing {fname}")

    available_stations = available_stations[
        available_stations["distance_km"] >= 0
    ].reset_index(drop=True)
    # + 10 because we expect some station with no data
    if len(available_stations) > (args.number_stations + 10):
        dmax = available_stations.iloc[args.number_stations + 10]["distance_km"]
        print(dmax)
        # Create a boolean mask and filter by distance
        mask = (available_stations["distance_km"] >= 0) & (
            available_stations["distance_km"] <= dmax
        )
        available_stations = available_stations[mask]

    print("available")
    print(available_stations)
    # required if not enough stations in the inventory
    args.number_stations = min(args.number_stations, len(available_stations))
    # initialize empty df
    selected_stations = pd.DataFrame(columns=available_stations.columns)


    while True:
        if args.closest_stations and selected_stations.empty:
            assert args.closest_stations <= args.number_stations
            previous_selected_stations = selected_stations.copy()
            selected_stations, available_stations = select_closest_stations(
                available_stations, selected_stations, args.closest_stations
            )
        else:
            if args.azimuthal:
                (
                    selected_stations,
                    available_stations,
                ) = select_telseismic_stations_aiming_for_azimuthal_coverage(
                    inventory, event, args.number_stations
                )
            else:
                previous_selected_stations = selected_stations.copy()
                selected_stations, available_stations = select_stations_most_distant(
                    available_stations, selected_stations, args.number_stations
                )

        added_rows = selected_stations[
            ~selected_stations["code"].isin(previous_selected_stations["code"])
        ]
        print("selection:", selected_stations)
        print("added:", added_rows)
        print("available:", available_stations)

        network_station = compute_dict_network_station(added_rows)

        retrieved_waveforms = retrieve_waveforms_including_preprocessed(
            network_station,
            client_name,
            kind_vd,
            path_observations,
            starttime,
            endtime,
            processed_data,
        )

        retrieved_stations = list(retrieved_waveforms.keys())
        print("retrieved_stations", retrieved_stations)
        added_rows = selected_stations[
            selected_stations["code"].isin(retrieved_stations)
        ]
        print("retrieved rows", added_rows)
        selected_stations = pd.concat(
            [previous_selected_stations, added_rows], ignore_index=True
        )
        print("new selected_stations", selected_stations)

        # required if not enough stations in the inventory
        args.number_stations = min(
            args.number_stations, len(selected_stations) + len(available_stations)
        )

        if len(selected_stations) == args.number_stations:
            print("done selecting stations")
            print(selected_stations)
            break
    generate_station_map(
        selected_stations, event, setup_name=setup_name, fault_info=fault_info
    )
    if "{{ stations }}" in config_stations:
        templateLoader = jinja2.FileSystemLoader(searchpath=os.getcwd())
        templateEnv = jinja2.Environment(loader=templateLoader)
        template = templateEnv.get_template(args.config_file)
        outputText = template.render({"stations": ",".join(selected_stations["code"])})
        out_fname = args.config_file
        with open(out_fname, "w") as fid:
            fid.write(outputText)
            print(f"done creating {out_fname}")
    else:
        print("no {{ stations }} field in the GENERAL/stations in the config file")
        print("please manually add:")
        station_codes = ",".join(selected_stations["code"])
        print(f"stations = {station_codes}")
