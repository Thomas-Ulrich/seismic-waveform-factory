#!/usr/bin/env python3
from obspy import UTCDateTime
from obspy.clients.fdsn import Client, RoutingClient
import pandas as pd
from obspy.geodetics.base import gps2dist_azimuth
import numpy as np
from plot_teleseismic_station_map import generate_station_map
import argparse
import configparser
import os
import jinja2
from retrieve_waveforms import retrieve_waveforms
from obspy import read_inventory
from geopy.distance import geodesic
from shapely.geometry import Point
import geopandas as gpd

np.random.seed(42)


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


def select_stations_most_distant(inv0, nstations):
    df = pd.DataFrame(
        columns=[
            "network",
            "station",
            "longitude",
            "latitude",
        ]
    )
    for i, net in enumerate(inv0):
        for j, sta in enumerate(net):
            new_row = {
                "network": net.code,
                "station": sta.code,
                "longitude": sta.longitude,
                "latitude": sta.latitude,
            }
            df.loc[len(df)] = new_row
    df["code"] = df["network"] + "." + df["station"]
    # Convert the latitude and longitude columns to a GeoDataFrame
    geometry = gpd.points_from_xy(df.longitude, df.latitude)
    crs = {"init": "epsg:4326"}
    gdf = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    gdf0 = gdf.copy()
    # Initialize with a random starting station
    selected_stations = gdf.sample(1)
    # Iterate through the remaining stations and select the one that is farthest from the currently selected stations
    for _ in range(nstations - 1):
        max_distance = 0
        furthest_station = None

        for index, row in gdf.iterrows():
            distance = min(
                [
                    haversine_distance(
                        row.latitude, row.longitude, station.latitude, station.longitude
                    )
                    for station in selected_stations.itertuples()
                ]
            )
            if distance > max_distance:
                max_distance = distance
                furthest_station = row
        selected_stations = pd.concat(
            [selected_stations, furthest_station.to_frame().T]
        )
        gdf = gdf.drop(
            index=furthest_station.name
        )  # Remove the selected station from the GeoDataFrame
    print(selected_stations)
    return selected_stations, gdf0


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


parser = argparse.ArgumentParser(
    description="select n stations ensuring an optimal coverage of an earthquake"
)
parser.add_argument("config_file", help="config file describing the event")
parser.add_argument("number_stations", help="number of stations to select", type=int)
parser.add_argument(
    "--azimuthal",
    dest="azimuthal",
    default=False,
    action="store_true",
    help="select stations based on back-azimuth (else based on distance)",
)

args = parser.parse_args()


config = configparser.ConfigParser()
assert os.path.isfile(args.config_file), f"{args.config_file} not found"
config.read(args.config_file)

client_name = config.get("GENERAL", "client")
onset = config.get("GENERAL", "onset")
hypo_lon = config.getfloat("GENERAL", "hypo_lon")
hypo_lat = config.getfloat("GENERAL", "hypo_lat")
stations = config.get("GENERAL", "stations")
path_observations = config.get("GENERAL", "path_observations")
kind_vd = config.get("GENERAL", "kind")
station_codes = stations.split(",")
software = config.get("GENERAL", "software").split(",")
is_teleseismic = "axitra" not in software

if not is_teleseismic:
    duration = config.getfloat("GENERAL", "axitra_duration")
    t_before = 100
    t_after = duration + 100
    min_max_radius = [1.0, 2.5]
else:
    import instaseis

    db_name = config.get("GENERAL", "db")
    db = instaseis.open_db(db_name)
    t_before = 1000
    t_after = db.info.length + 1000
    min_max_radius = [30, 90]


if "{{ stations }}" in stations:
    if client_name in ["eida-routing", "iris-federator"]:
        client = RoutingClient(client_name)
    else:
        client = Client(client_name)

    # Define the event parameters
    event_time = UTCDateTime(onset)
    event = {"longitude": hypo_lon, "latitude": hypo_lat, "onset": event_time}
    starttime = event["onset"] - t_before
    endtime = event["onset"] + t_after

    fn_inventory = "observations/inv_{client_name}.xml"
    if os.path.exists(fn_inventory):
        inventory = read_inventory(fn_inventory)
    else:
        # Get station information around the event
        inventory = client.get_stations(
            latitude=event["latitude"],
            longitude=event["longitude"],
            starttime=starttime,
            endtime=endtime,
            level="channel",
            # includerestricted=False,
            # includeavailability = True,
            # matchtimeseries=True,
            minradius=min_max_radius[0],
            maxradius=min_max_radius[1],
            channel="*Z",
        )
        inventory.write(fn_inventory, format="STATIONXML")
    print(inventory)

    if not args.azimuthal:
        df, df_all = select_stations_most_distant(inventory, args.number_stations)
    else:
        df, df_all = select_telseismic_stations_aiming_for_azimuthal_coverage(
            inventory, event, args.number_stations
        )

    network_station = compute_dict_network_station(df)
    retrieved_waveforms = retrieve_waveforms(
        network_station, client_name, kind_vd, path_observations, starttime, endtime
    )
    retrieved_stations = list(retrieved_waveforms.keys())
    index_valid_stations = []
    while len(retrieved_stations) < args.number_stations:
        index_new_stations = []
        for index, row in df.iterrows():
            row_index = df_all.index[df_all["code"] == row["code"]].tolist()
            if not row["code"] in retrieved_stations:
                index_new_stations.append(row_index[0] + 1)
            else:
                index_valid_stations.append(row_index[0])
        df = df_all.loc[index_new_stations]
        print(df)
        network_station = compute_dict_network_station(df)
        new_retrieved_waveforms = retrieve_waveforms(
            network_station, client_name, kind_vd, path_observations, starttime, endtime
        )
        retrieved_stations.extend(new_retrieved_waveforms.keys())

        print(retrieved_stations)
    # regenerate the df of all the valid stations
    for index, row in df.iterrows():
        row_index = df_all.index[df_all["code"] == row["code"]].tolist()
        index_valid_stations.append(row_index[0])

    df = df_all.loc[index_valid_stations]
    print(df)

    generate_station_map(df, event)

    templateLoader = jinja2.FileSystemLoader(searchpath=os.getcwd())
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(args.config_file)
    outputText = template.render({"stations": ",".join(df["code"])})
    out_fname = args.config_file
    with open(out_fname, "w") as fid:
        fid.write(outputText)
        print(f"done creating {out_fname}")

else:
    print(f"station_codes is already set to {stations} in {args.config_file}")
