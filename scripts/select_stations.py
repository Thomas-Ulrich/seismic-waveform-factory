#!/usr/bin/env python3
import argparse
import configparser
import os

import geopandas as gpd
import jinja2
import numpy as np
import pandas as pd
from fault_processing import compute_shapely_polygon, get_fault_slip_coords
from geodetic_utils import add_distance_backazimuth_to_df
from geopy.distance import geodesic
from obspy import UTCDateTime, read_inventory
from obspy.clients.fdsn import Client, RoutingClient
from obspy.core.inventory import Inventory
from plot_station_map import generate_station_map
from pyproj import Transformer
from retrieve_waveforms import filter_channels_by_availability, retrieve_waveforms
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


# Define a function to calculate the distance between
# two points using their latitude and longitude
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


def generate_geopanda_dataframe(df_stations, fault_info, projection):
    df = pd.DataFrame(
        columns=["network", "station", "longitude", "latitude", "distance_km"]
    )
    if projection:
        transformer = Transformer.from_crs("epsg:4326", projection, always_xy=True)
        if fault_info:
            coords = fault_info["fault_slip_coords"] / 1e3
            tree = spatial.KDTree(coords)
        else:
            print("using hypocenter to compute distance in km")
            x1, y1 = transformer.transform(hypo_lon, hypo_lat)
            hypo_coords = np.array([x1 / 1e3, y1 / 1e3, -hypo_depth_in_km])

    for _, row in df_stations.iterrows():
        if projection:
            x1, y1 = transformer.transform(row.longitude, row.latitude)
            if fault_info:
                min_distance, _ = tree.query([x1 / 1e3, y1 / 1e3, 0])
            else:
                min_distance = np.linalg.norm(
                    np.array([x1 / 1e3, y1 / 1e3, 0]) - hypo_coords
                )
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


def select_teleseismic_stations_aiming_for_azimuthal_coverage(
    available_stations, selected_stations_at_start, nstations_to_select
):
    nstations = nstations_to_select - len(selected_stations_at_start)
    df = available_stations
    # Parameters
    n = 8  # Number of backazimuth panels (0-360 degrees)
    p = 4  # Number of distance panels (0-90 degrees)

    # Function to assign backazimuth panel
    def assign_backazimuth_panel(backazimuth, n):
        return int(backazimuth // (360 / n))

    # Function to assign distance panel
    def assign_distance_panel(distance, p):
        return int((distance - 30) // (60 / p))

    # Assign panels based on backazimuth and distance
    df["backazimuth_panel"] = df["backazimuth"].apply(assign_backazimuth_panel, n=n)
    df["distance_panel"] = df["distance"].apply(assign_distance_panel, p=p)

    # Group stations by their backazimuth and distance panels
    panel_groups = df.groupby(["backazimuth_panel", "distance_panel"])

    # Total number of panels available
    num_panels = len(panel_groups)

    selected_stations = []

    if num_panels >= nstations:
        # If there are enough panels, select `n` random panels and sample one station
        # per panel
        # Convert panel_groups.groups.keys() to a list
        panel_keys = list(panel_groups.groups.keys())

        # Randomly select panels
        selected_panels = np.random.choice(
            len(panel_keys), size=nstations, replace=False
        )

        # Map the selected indices back to the actual keys
        selected_panels = [panel_keys[i] for i in selected_panels]

        for panel in selected_panels:
            group = panel_groups.get_group(panel)
            selected_station = group.sample(n=1)  # Randomly select one station
            selected_stations.append(selected_station)
        selected_stations_df = pd.concat(selected_stations)

    else:
        # If there are fewer panels than `n`, select all available panels
        selected_stations = []
        for (backazimuth_panel, distance_panel), group in panel_groups:
            # Randomly select one station from each group
            selected_station = group.sample(n=1)
            selected_stations.append(selected_station)

        # Combine the selected stations into a DataFrame to track them
        selected_stations_df = pd.concat(selected_stations)
        selected_station_ids = selected_stations_df["station"].values.tolist()

        # If there are still remaining stations to reach `nstations`,
        # sample from the remaining stations
        remaining_stations_needed = nstations - len(selected_stations_df)
        if remaining_stations_needed > 0:
            remaining_groups = panel_groups.apply(
                lambda x: x[~x["station"].isin(selected_station_ids)]
            )
            remaining_stations = remaining_groups.sample(
                n=remaining_stations_needed, replace=False
            )

            selected_stations_df = pd.concat([selected_stations_df, remaining_stations])

    # Final DataFrame of selected stations
    df_selected = pd.merge(
        selected_stations_df, selected_stations_at_start, how="outer"
    )

    df_remaining = df[~df["station"].isin(df_selected["station"])]
    return df_selected, df_remaining


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
        "--distance_range",
        type=float,
        help="distance range from selecting stations (degrees)",
        nargs=2,
    )

    parser.add_argument(
        "--channel", default="*", type=str, help="filter channels to be retrieved"
    )
    parser.add_argument(
        "--store_format",
        choices=["sac", "mseed"],
        default="mseed",
        type=str,
        help="store format for waveform data",
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
    client,
    client_name,
    event,
    path_observations,
    spatial_range,
    t_before,
    t_after,
    channel,
):
    fn_inventory = f"{path_observations}/inv_{client_name}.xml"
    if os.path.exists(fn_inventory):
        inventory = read_inventory(fn_inventory)
    else:
        if channel != "*":
            print(f"filtering channels, channel = {channel}")
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
            if r1 > 30:
                kargs["network"] = "IU,II,GE,G"
                print(
                    f"Warning: using only networks {kargs['network']} for teleseismic"
                )

        if client_name in ["NCEDC"]:
            level = "station"
            print(
                f"{client_name} won't return channel information for multiple station\
            using level = {level}"
            )
        else:
            level = "channel"
            kargs["includerestricted"] = False

        inventory = client.get_stations(
            starttime=starttime,
            endtime=endtime,
            level=level,
            channel=channel,
            includeavailability=True,
            **kargs,
        )

        if level == "channel":
            inventory = filter_channels_by_availability(inventory, starttime, endtime)
        inventory.write(fn_inventory, format="STATIONXML")

    inventory = inventory.select(channel=channel)
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
    hypo_depth_in_km = config.getfloat("GENERAL", "hypo_depth_in_km")
    config_stations = config.get("GENERAL", "stations", fallback="")
    path_observations = config.get("GENERAL", "path_observations")
    kind_vd = config.get("GENERAL", "kind")
    software = config.get("GENERAL", "software").split(",")
    is_teleseismic = "axitra" not in software

    faultfname = config.get("GENERAL", "fault_filename", fallback=None)
    station_file = config.get("GENERAL", "station_file", fallback=None)

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
        try:
            import instaseis

            db_name = config.get("GENERAL", "db")
            db = instaseis.open_db(db_name)
            signal_length = db.info.length
        except (ImportError, ModuleNotFoundError) as e:
            print(f"instaseis not available: {e}")
            signal_length = 3600.0
        except (OSError, IOError, ValueError) as e:
            print(f"Could not open Instaseis DB: {e}")
            signal_length = 3600.0

        t_before = 1000
        t_after = signal_length + 1000
        spatial_range["radius"] = [30.0, 90.0]

    if args.distance_range:
        spatial_range["radius"] = args.distance_range

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
    else:
        inventory = load_or_create_inventory(
            client,
            client_name,
            event,
            path_observations,
            spatial_range,
            t_before,
            t_before,
            args.channel,
        )
        inventory = remove_synthetics_from_inventory(inventory)
        filtered_networks = [net for net in inventory.networks if net.code != "AM"]
        new_inventory = Inventory(networks=filtered_networks, source=inventory.source)
        all_stations = new_inventory.get_contents()["stations"]
        if len(all_stations) > (args.number_stations):
            print("remove AM network")
            inventory = new_inventory
        else:
            print("did not remove AM network, else too no enough stations remaining")

        station_df = generate_station_df(inventory)
    print(station_df)

    projection = config.get("GENERAL", "projection", fallback="")
    available_stations = generate_geopanda_dataframe(station_df, fault_info, projection)
    print(available_stations)

    if projection:
        # 60 closest stations are written for seissol output
        closest_stations = available_stations.head(60)
        transformer = Transformer.from_crs("epsg:4326", projection, always_xy=True)
        x1, y1 = transformer.transform(
            closest_stations["longitude"], closest_stations["latitude"]
        )
        n_seissol_station = len(closest_stations)
        nstations_in_file = [n_seissol_station]
        if n_seissol_station > 5:
            nstations_in_file.append(5)

        for k in nstations_in_file:
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
    if len(available_stations) > (args.number_stations + 10) and not is_teleseismic:
        dmax = available_stations.iloc[args.number_stations + 10]["distance_km"]
        # Create a boolean mask and filter by distance
        mask = (available_stations["distance_km"] >= 0) & (
            available_stations["distance_km"] <= dmax
        )
        other_available_stations = available_stations[~mask]
        available_stations = available_stations[mask]
        print(
            f"available ({args.number_stations + 10} closest, that is up to {dmax} km):"
        )
    else:
        other_available_stations = gpd.GeoDataFrame()
        print("available (no restrictions):")

    if args.azimuthal:
        available_stations = add_distance_backazimuth_to_df(available_stations, event)
        available_stations = available_stations.drop(columns=["geometry"])

    print(available_stations)
    # required if not enough stations in the inventory
    args.number_stations = min(args.number_stations, len(available_stations))
    # initialize empty df
    selected_stations = pd.DataFrame(columns=available_stations.columns)

    while True:
        if args.closest_stations and selected_stations.empty:
            if args.closest_stations <= args.number_stations:
                print("args.closest_stations <= args.number_stations")
                args.closest_stations = args.number_stations
            previous_selected_stations = selected_stations.copy()
            selected_stations, available_stations = select_closest_stations(
                available_stations, selected_stations, args.closest_stations
            )
        else:
            previous_selected_stations = selected_stations.copy()
            if args.azimuthal:
                (
                    selected_stations,
                    available_stations,
                ) = select_teleseismic_stations_aiming_for_azimuthal_coverage(
                    available_stations, selected_stations, args.number_stations
                )
            else:
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
        # transform the dictionnary in a list of strings "network.station"
        network_station = [
            f"{key}.{value}"
            for key, values in network_station.items()
            for value in values
        ]

        retrieved_waveforms = retrieve_waveforms(
            network_station,
            client_name,
            kind_vd,
            path_observations,
            starttime,
            endtime,
            processed_data=processed_data,
            output_format=args.store_format,
            channel=args.channel,
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
        if not other_available_stations.empty and len(
            available_stations
        ) < args.number_stations - len(selected_stations):
            available_stations = gpd.GeoDataFrame(
                pd.concat(
                    [available_stations, other_available_stations], ignore_index=True
                )
            )
            print("print adding first discarded other available stations")
            print(other_available_stations)
            print("available_stations is now:")
            print(available_stations)
            other_available_stations = gpd.GeoDataFrame()

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
