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


def generate_geopanda_dataframe(inv0):
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
    crs = "epsg:4326"
    return gpd.GeoDataFrame(df, crs=crs, geometry=geometry)


def select_stations_most_distant(available_stations, selected_stations, nstations):
    if len(selected_stations) == 0:
        print("selected_stations is empty, selecting randomly from available_stations")
        # Initialize with a random starting station
        # and remove the selected station from the available_stations
        selected_stations = available_stations.sample(1)
        available_stations = available_stations.drop(
            index=selected_stations.index.values[0]
        )

    # Iterate through the remaining stations and select the one that is farthest from the currently selected stations
    nselected = len(selected_stations)
    for k in range(nstations - nselected):
        if available_stations.empty:
            print(f"GeoDataFrame is empty, {nstations-1-k} stations will be missing")
            break

        max_distance = 0
        furthest_station = None

        for index, row in available_stations.iterrows():
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
        available_stations = available_stations.drop(
            index=furthest_station.name
        )  # Remove the selected station from the GeoDataFrame
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="select n stations ensuring an optimal coverage of an earthquake"
    )
    parser.add_argument("config_file", help="config file describing the event")
    parser.add_argument(
        "number_stations", help="number of stations to select", type=int
    )
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

    if config.has_option("GENERAL", "stations"):
        config_stations = config.get("GENERAL", "stations")
    else:
        config_stations = ""

    path_observations = config.get("GENERAL", "path_observations")
    kind_vd = config.get("GENERAL", "kind")
    software = config.get("GENERAL", "software").split(",")
    is_teleseismic = "axitra" not in software

    if not is_teleseismic:
        duration = config.getfloat("GENERAL", "axitra_pyprop8_duration")
        t_before = 100
        t_after = duration + 100
        min_max_radius = [0.0, 2.5]
    else:
        import instaseis

        db_name = config.get("GENERAL", "db")
        db = instaseis.open_db(db_name)
        t_before = 1000
        t_after = db.info.length + 1000
        min_max_radius = [30, 90]

    if client_name in ["eida-routing", "iris-federator"]:
        client = RoutingClient(client_name)
    else:
        client = Client(client_name)

    # Define the event parameters
    event_time = UTCDateTime(onset)
    event = {"longitude": hypo_lon, "latitude": hypo_lat, "onset": event_time}
    starttime = event["onset"] - t_before
    endtime = event["onset"] + t_after

    if not os.path.exists(path_observations):
        os.makedirs(path_observations)
        print(f"done creating {path_observations}")

    fn_inventory = f"{path_observations}/inv_{client_name}.xml"
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

    available_stations = generate_geopanda_dataframe(inventory)
    print('available')
    print(available_stations)
    # required if not enough stations in the inventory
    args.number_stations = min(args.number_stations, len(available_stations))
    # initialize empty df
    selected_stations = pd.DataFrame(columns=available_stations.columns)

    while True:

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

        added_rows = selected_stations[~selected_stations['code'].isin(previous_selected_stations['code'])]
        print("selection:", selected_stations)
        print("added:", added_rows)
        print("available:", available_stations)

        network_station = compute_dict_network_station(added_rows)
        retrieved_waveforms = retrieve_waveforms(
            network_station, client_name, kind_vd, path_observations, starttime, endtime
        )
        retrieved_stations = list(retrieved_waveforms.keys())
        print("retrieved_stations", retrieved_stations)
        added_rows = selected_stations[selected_stations['code'].isin(retrieved_stations)]
        print("retrieved rows", added_rows)
        selected_stations = pd.concat([previous_selected_stations, added_rows], ignore_index=True)
        print("new selected_stations", selected_stations)
        
        # required if not enough stations in the inventory
        args.number_stations = min(args.number_stations, len(selected_stations) + len(available_stations))

        if len(selected_stations) == args.number_stations:
            print("done selecting stations")
            print(selected_stations)
            break

    generate_station_map(selected_stations, event)
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
