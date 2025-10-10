#!/usr/bin/env python3
import argparse
import os
import random

import geopandas as gpd
import jinja2
import numpy as np
import pandas as pd
from config_loader import ConfigLoader
from config_schema import CONFIG_SCHEMA
from config_utils import (
    categorize_stations_by_scale,
    categorize_waveform_kind_by_scale,
    extract_instaseis_db,
    extract_regional_durations,
)
from fault_processing import compute_shapely_polygon, get_fault_slip_coords
from geodetic_utils import add_distance_backazimuth_to_df
from geopy.distance import geodesic
from obspy import UTCDateTime, read_inventory
from obspy.core.inventory import Inventory
from plot_station_map import generate_station_map
from pyproj import Transformer
from retrieve_waveforms import (
    filter_channels_by_availability,
    initialize_client,
    retrieve_waveforms,
)
from scipy import spatial

# from waveform_figure_utils import get_station_files_dict

np.random.seed(42)
random.seed(42)


def generate_station_df(inv0):
    data = [
        {
            "network": net.code,
            "station": sta.code,
            "longitude": sta.longitude,
            "latitude": sta.latitude,
        }
        for net in inv0
        for sta in net
    ]
    return pd.DataFrame(data)


def compute_dict_network_station(df):
    networks = set(df["network"])
    network_station = {}
    for net in networks:
        network_station[net] = []
    for index, row in df.iterrows():
        network_station[row["network"]].append(row["station"])
    return network_station


def haversine_distance(lat1, lon1, lat2, lon2):
    """Return the geodesic distance (km) between two latitude-longitude points."""
    return geodesic((lat1, lon1), (lat2, lon2)).km


def remove_synthetics_from_inventory(original_inv):
    """Return a copy of the inventory without synthetic stations."""
    new_inv = Inventory()
    for net in original_inv:
        stations = [
            sta
            for sta in net
            if not sta.site.name or "synthetic" not in sta.site.name.lower()
        ]
        if stations:
            new_net = net.copy()
            new_net.stations = stations
            new_inv.networks.append(new_net)
    return new_inv


def generate_geopanda_dataframe(df_stations, event, fault_info=None, projection=None):
    """Return a GeoDataFrame of stations with distances (km) to fault or hypocenter."""
    records = []
    print(df_stations, event, fault_info, projection)
    transformer = (
        Transformer.from_crs("epsg:4326", projection, always_xy=True)
        if projection
        else None
    )

    # Setup distance reference (fault or hypocenter)
    if projection and fault_info:
        tree = spatial.KDTree(fault_info["fault_slip_coords"] / 1e3)
    elif projection:
        print("Using hypocenter to compute distance in km")
        x1, y1 = transformer.transform(event["lon"], event["lat"])
        hypo_coords = np.array([x1, y1, -event["hypo_depth_in_km"]]) / 1e3
    else:
        tree = hypo_coords = None

    # Compute distances
    for _, row in df_stations.iterrows():
        if not projection:
            dist = -1.0
        else:
            x, y = transformer.transform(row.longitude, row.latitude)
            print(x, y)
            xyz = np.array([x, y, 0]) / 1e3
            dist = (
                tree.query(xyz)[0] if fault_info else np.linalg.norm(xyz - hypo_coords)
            )

        records.append(
            {
                "network": row.network,
                "station": row.station,
                "longitude": row.longitude,
                "latitude": row.latitude,
                "distance_km": dist,
            }
        )

    df = pd.DataFrame(records)
    df["code"] = df["network"] + "." + df["station"]
    df.sort_values("distance_km", inplace=True, ignore_index=True)

    return gpd.GeoDataFrame(
        df,
        crs="epsg:4326",
        geometry=gpd.points_from_xy(df.longitude, df.latitude),
    )


def select_closest_stations(available_stations, selected_stations, nstations):
    """Return updated DataFrames of selected and remaining stations."""
    if selected_stations.empty:
        print("selected_stations is empty; selecting the closest stations.")
        selected_stations = available_stations.head(nstations)
        available_stations = available_stations.iloc[nstations:]

    return selected_stations.sort_index(), available_stations.sort_index()


def select_stations_maximizing_distance(
    available_stations, selected_stations, nstations
):
    """Select stations maximizing distance diversity until reaching nstations."""
    if selected_stations.empty:
        print("selected_stations is empty; selecting one station at random.")
        selected_stations = available_stations.sample(1)
        available_stations = available_stations.drop(selected_stations.index)

    while len(selected_stations) < nstations and not available_stations.empty:
        # Compute minimum distance of each available station to any selected ones
        distances = available_stations.apply(
            lambda row: min(
                haversine_distance(row.latitude, row.longitude, s.latitude, s.longitude)
                for s in selected_stations.itertuples(index=False)
            ),
            axis=1,
        )

        # Pick the farthest station and update sets
        furthest_idx = distances.idxmax()
        selected_stations = pd.concat(
            [selected_stations, available_stations.loc[[furthest_idx]]]
        )
        available_stations = available_stations.drop(furthest_idx)

    if available_stations.empty and len(selected_stations) < nstations:
        print("GeoDataFrame is empty;")
        print(f"{nstations - len(selected_stations)} stations missing.")

    return selected_stations.sort_index(), available_stations.sort_index()


def select_teleseismic_stations_aiming_for_azimuthal_coverage(
    available_stations, selected_stations_at_start, nstations_to_select
):
    """
    Select teleseismic stations aiming for even azimuthal and distance coverage.

    Parameters
    ----------
    available_stations : pd.DataFrame
        Stations with 'backazimuth' and 'distance' columns.
    selected_stations_at_start : pd.DataFrame
        Already selected stations to include.
    nstations_to_select : int
        Total number of stations to select.

    Returns
    -------
    pd.DataFrame, pd.DataFrame
        Selected stations, remaining stations.
    """
    nstations = nstations_to_select - len(selected_stations_at_start)
    if nstations <= 0:
        return selected_stations_at_start.copy(), available_stations.copy()

    df = available_stations.copy()

    # Define panels
    n_baz_panels = 8
    n_dist_panels = 4

    # Assign panels
    df["backazimuth_panel"] = (df["backazimuth"] // (360 / n_baz_panels)).astype(int)
    df["distance_panel"] = ((df["distance"] - 30) // (60 / n_dist_panels)).astype(int)

    panel_groups = df.groupby(["backazimuth_panel", "distance_panel"])
    panel_keys = list(panel_groups.groups.keys())

    selected_stations_list = []

    if len(panel_keys) >= nstations:
        chosen_panels = random.sample(panel_keys, k=nstations)
        for panel in chosen_panels:
            selected_stations_list.append(panel_groups.get_group(panel).sample(n=1))
    else:
        # Select one station per panel
        for _, group in panel_groups:
            selected_stations_list.append(group.sample(n=1))

        selected_ids = pd.concat(selected_stations_list)["station"].tolist()
        remaining_needed = nstations - len(selected_ids)

        if remaining_needed > 0:
            remaining_stations = df[~df["station"].isin(selected_ids)]
            selected_stations_list.append(
                remaining_stations.sample(n=remaining_needed, replace=False)
            )

    selected_stations_df = pd.concat(selected_stations_list)
    df_selected = pd.merge(
        selected_stations_at_start, selected_stations_df, how="outer"
    )
    df_remaining = df[~df["station"].isin(df_selected["station"])]

    return df_selected.reset_index(drop=True), df_remaining.reset_index(drop=True)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Select stations ensuring optimal coverage for an earthquake.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional arguments
    parser.add_argument("config_file", help="yaml config file describing the event.")
    parser.add_argument(
        "number_stations", type=int, help="Total number of stations to select."
    )
    parser.add_argument(
        "closest_stations", type=int, help="Number of closest stations to select."
    )

    # Optional arguments
    parser.add_argument(
        "--distance_range",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="Distance range (degrees) from which to select stations.",
    )
    parser.add_argument(
        "--channel",
        default="*",
        type=str,
        help='Filter channels to be retrieved (default "*" = all channels).',
    )
    parser.add_argument(
        "--store_format",
        choices=["sac", "mseed"],
        default="mseed",
        type=str,
        help="Storage format for waveform data.",
    )
    parser.add_argument(
        "--azimuthal",
        action="store_true",
        help="Select stations based on back-azimuth instead of distance.",
    )

    return parser.parse_args()


def normalize_bounds(spatial_range, r1):
    """Expand spatial bounds and normalize lat/lon within valid Earth ranges."""

    def clamp_lon(lon):
        """Clamp longitude to [-180, 180]."""
        return max(-180, min(180, lon))

    def clamp_lat(lat):
        """Clamp latitude to [-90, 90]."""
        return max(-90, min(90, lat))

    lon_min = clamp_lon(spatial_range["source_lon_range"][0] - r1)
    lon_max = clamp_lon(spatial_range["source_lon_range"][1] + r1)
    lat_min = clamp_lat(spatial_range["source_lat_range"][0] - r1)
    lat_max = clamp_lat(spatial_range["source_lat_range"][1] + r1)

    kargs = {
        "minlongitude": lon_min,
        "maxlongitude": lon_max,
        "minlatitude": lat_min,
        "maxlatitude": lat_max,
    }

    return kargs


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
        if "source_lon_range" in spatial_range.keys() and r1 < 30:
            kargs = normalize_bounds(spatial_range, r1)
        else:
            kargs["minradius"] = r0
            kargs["maxradius"] = r1
            kargs["latitude"] = event["lat"]
            kargs["longitude"] = event["lon"]
            if r1 >= 30:
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
        kargs["starttime"] = starttime
        kargs["endtime"] = endtime
        kargs["level"] = level
        kargs["channel"] = channel
        kargs["includeavailability"] = True

        print("running client.get_stations with:")
        print(kargs)
        inventory = client.get_stations(
            **kargs,
        )
        print("done running client.get_stations")

        if level == "channel":
            inventory = filter_channels_by_availability(inventory, starttime, endtime)
        inventory.write(fn_inventory, format="STATIONXML")

    inventory = inventory.select(channel=channel)
    return inventory


def compute_min_max_coords(fault_info):
    transformer = Transformer.from_crs(
        fault_info["projection"], "epsg:4326", always_xy=True
    )
    xs, ys = [], []

    for poly in fault_info["polygons"]:
        x, y = transformer.transform(*poly.exterior.xy)
        xs.extend(x)
        ys.extend(y)

    return (min(xs), max(xs)), (min(ys), max(ys))


def select_station(
    config_file,
    number_stations,
    closest_stations,
    distance_range,
    channel,
    store_format,
    azimuthal,
):
    cfg = ConfigLoader(config_file, CONFIG_SCHEMA)

    client_name = cfg["general"]["client"]
    path_observations = cfg["general"]["path_observations"]

    station_codes = categorize_stations_by_scale(cfg)
    waveform_kind = categorize_waveform_kind_by_scale(cfg)

    is_teleseismic = True if station_codes["global"] else False
    if station_codes["global"] and station_codes["regional"]:
        print(
            "there are plots with and without instaseis.",
            "Selecting only teleseismic stations",
        )

    key = "global" if is_teleseismic else "regional"
    kind_vd = waveform_kind[key][0]
    assert (
        len(waveform_kind[key]) == 1
    ), f"several waveform kind requested, {waveform_kind[key]}"

    station_file = cfg["general"]["station_file"]

    processed_data = {}
    """
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

    spatial_range = {}

    faultfname = cfg["general"]["fault_file"]
    if faultfname:
        polygons = compute_shapely_polygon(faultfname)
        fault_slip_coords = get_fault_slip_coords(faultfname)
        projection = cfg["general"]["projection"]
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
        durations = extract_regional_durations(cfg)
        signal_length = max(durations)
        spatial_range["radius"] = [0.0, 2.5]
    else:
        try:
            import instaseis

            dbs = extract_instaseis_db(cfg)
            signal_length = 0
            for db_name in dbs:
                db = instaseis.open_db(db_name)
                signal_length = max(signal_length, db.info.length)
        except (ImportError, ModuleNotFoundError) as e:
            print(f"instaseis not available: {e}")
            signal_length = 3600.0
        except (OSError, IOError, ValueError) as e:
            print(f"Could not open Instaseis DB: {e}")
            signal_length = 3600.0

        spatial_range["radius"] = [30.0, 90.0]

    extra_time = max(100.0, 0.1 * signal_length)
    t_before = extra_time
    t_after = signal_length + extra_time

    if distance_range:
        spatial_range["radius"] = distance_range

    client = initialize_client(client_name)

    # Define the event parameters
    event = cfg["general"]["hypocenter"]
    event["onset"] = UTCDateTime(event["onset"])

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
            channel,
        )
        inventory = remove_synthetics_from_inventory(inventory)
        filtered_networks = [net for net in inventory.networks if net.code != "AM"]
        new_inventory = Inventory(networks=filtered_networks, source=inventory.source)
        all_stations = new_inventory.get_contents()["stations"]
        if len(all_stations) > (number_stations):
            print("remove AM network")
            inventory = new_inventory
        else:
            print("did not remove AM network, else too no enough stations remaining")

        station_df = generate_station_df(inventory)
    print(station_df)

    projection = cfg["general"]["projection"]
    if is_teleseismic:
        projection = None

    available_stations = generate_geopanda_dataframe(
        station_df, event, fault_info, projection
    )
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
    if len(available_stations) > (number_stations + 10) and not is_teleseismic:
        dmax = available_stations.iloc[number_stations + 10]["distance_km"]
        # Create a boolean mask and filter by distance
        mask = (available_stations["distance_km"] >= 0) & (
            available_stations["distance_km"] <= dmax
        )
        other_available_stations = available_stations[~mask]
        available_stations = available_stations[mask]
        print(f"available ({number_stations + 10} closest, that is up to {dmax} km):")
    else:
        other_available_stations = gpd.GeoDataFrame()
        print("available (no restrictions):")

    if azimuthal:
        available_stations = add_distance_backazimuth_to_df(available_stations, event)
        available_stations = available_stations.drop(columns=["geometry"])

    print(available_stations)
    # required if not enough stations in the inventory
    number_stations = min(number_stations, len(available_stations))
    # initialize empty df
    selected_stations = pd.DataFrame(columns=available_stations.columns)

    while True:
        if closest_stations and selected_stations.empty:
            if closest_stations <= number_stations:
                print("closest_stations <= number_stations")
                closest_stations = number_stations
            previous_selected_stations = selected_stations.copy()
            selected_stations, available_stations = select_closest_stations(
                available_stations, selected_stations, closest_stations
            )
        else:
            previous_selected_stations = selected_stations.copy()
            if azimuthal:
                (
                    selected_stations,
                    available_stations,
                ) = select_teleseismic_stations_aiming_for_azimuthal_coverage(
                    available_stations, selected_stations, number_stations
                )
            else:
                (
                    selected_stations,
                    available_stations,
                ) = select_stations_maximizing_distance(
                    available_stations, selected_stations, number_stations
                )

        added_rows = selected_stations[
            ~selected_stations["code"].isin(previous_selected_stations["code"])
        ]
        print("selection:\n", selected_stations)
        print("added:\n", added_rows)
        print("available:\n", available_stations)

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
            output_format=store_format,
            channel=channel,
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
        ) < number_stations - len(selected_stations):
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

        number_stations = min(
            number_stations, len(selected_stations) + len(available_stations)
        )

        if len(selected_stations) == number_stations:
            print("done selecting stations")
            print(selected_stations)
            break

    generate_station_map(selected_stations, cfg, is_teleseismic)

    if "{{ stations }}" in config_stations:
        templateLoader = jinja2.FileSystemLoader(searchpath=os.getcwd())
        templateEnv = jinja2.Environment(loader=templateLoader)
        template = templateEnv.get_template(config_file)
        outputText = template.render({"stations": ",".join(selected_stations["code"])})
        out_fname = config_file
        with open(out_fname, "w") as fid:
            fid.write(outputText)
            print(f"done creating {out_fname}")
    else:
        print("no {{ stations }} field in the GENERAL/stations in the config file")
        print("please manually add:")
        station_codes = ",".join(selected_stations["code"])
        print(f"stations = {station_codes}")


if __name__ == "__main__":
    args = parse_arguments()
    select_station(
        args.config_file,
        args.number_stations,
        args.closest_stations,
        args.distance_range,
        args.channel,
        args.store_format,
        args.azimuthal,
    )
