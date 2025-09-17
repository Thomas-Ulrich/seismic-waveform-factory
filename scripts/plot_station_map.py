#!/usr/bin/env python3
from obspy import UTCDateTime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse
import configparser
import os
from fault_processing import compute_shapely_polygon
from retrieve_waveforms import initialize_client
from obspy import read_inventory
from obspy.core.inventory import Inventory
from geodetic_utils import add_distance_backazimuth_to_df
from scalebar import scale_bar
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import glob


def retrieve_coordinates(client_name, event, station_codes):
    df = pd.DataFrame(columns=["network", "station", "longitude", "latitude"])
    fn_inventories = glob.glob(os.path.join(path_observations, "inv_*.xml"))
    for fn_inventory in fn_inventories:
        inventory = read_inventory(fn_inventory)
        for network in inventory:
            for station in network:
                code = f"{network.code}.{station.code}"
                if code in station_codes:
                    new_row = {
                        "network": network.code,
                        "station": station.code,
                        "longitude": station.longitude,
                        "latitude": station.latitude,
                    }
                    df.loc[len(df)] = new_row

    else:
        client = initialize_client(client_name)

        event_time = UTCDateTime(event["onset"])
        networks = set([entry.split(".")[0] for entry in station_codes])
        print(networks)
        inv = Inventory()
        for net in networks:
            fn_inventory = f"{path_observations}/inv_stations_{net}.xml"
            if os.path.exists(fn_inventory):
                inventory = read_inventory(fn_inventory)
            else:
                inventory = client.get_stations(
                    network=net,
                    level="station",
                    starttime=event_time,
                    endtime=event_time + 100,
                    includeavailability=True,
                )
                inventory.write(fn_inventory, format="STATIONXML")
            inv.extend(inventory)
        network_station_pairs = [code.split(".") for code in station_codes]
        extracted_data = [
            (network.code, station.code, station.longitude, station.latitude)
            for network in inv.networks
            for station in network.stations
            if [network.code, station.code] in network_station_pairs
        ]
        print(extracted_data)
        for vals in extracted_data:
            net, sta, lon, lat = vals
            new_row = {
                "network": net,
                "station": sta,
                "longitude": lon,
                "latitude": lat,
            }
            df.loc[len(df)] = new_row

    df = add_distance_backazimuth_to_df(df, event)
    return df


def generate_station_map(df, event, set_global=False, setup_name="", fault_info=None):
    plt.figure(figsize=(6, 6))
    if set_global:
        projection = ccrs.Orthographic(
            central_longitude=event["longitude"], central_latitude=event["latitude"]
        )
        geo = ccrs.Geodetic()
        ax = plt.axes(projection=projection)
        scale = "110m"
    else:
        geo = ccrs.Geodetic()
        ax = plt.axes(projection=ccrs.PlateCarree())
        scale = "10m"

    ax.add_feature(cfeature.LAND.with_scale(scale), rasterized=True)
    ax.add_feature(cfeature.OCEAN.with_scale(scale), rasterized=True)
    ax.add_feature(cfeature.COASTLINE.with_scale(scale))
    ax.add_feature(cfeature.BORDERS.with_scale(scale), linestyle=":")
    x = df["longitude"].values
    y = df["latitude"].values

    if set_global:
        x1, x2 = min(event["longitude"], x.min()), max(event["longitude"], x.max())
        y1, y2 = min(event["latitude"], y.min()), max(event["latitude"], y.max())
        dx = 0.2 * (x2 - x1)
        dy = 0.2 * (y2 - y1)
        ax.set_extent(
            [
                max(-180, x1 - dx),
                min(180, x2 + dx),
                max(-90, y1 - dy),
                min(90, y2 + dy),
            ],
            ccrs.PlateCarree(),
        )
    else:
        # Add Latitude and Longitude Gridlines
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = mticker.FixedLocator(
            np.arange(-180, 181, 0.5)
        )  # Adjust as needed
        gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 0.5))  # Adjust as needed

    names = df["station"].values
    plt.scatter(
        x, y, 200, color="r", marker="v", edgecolor="k", zorder=3, transform=geo
    )
    plt.scatter(
        event["longitude"],
        event["latitude"],
        200,
        color="orange",
        marker="*",
        edgecolor="k",
        zorder=3,
        transform=geo,
    )
    if fault_info:
        from pyproj import Transformer

        projection = fault_info["projection"]
        transformer = Transformer.from_crs(projection, "epsg:4326", always_xy=True)
        polygons = fault_info["polygons"]
        for poly in polygons:
            x1, y1 = poly.exterior.xy
            x1, y1 = transformer.transform(x1, y1)
            plt.plot(x1, y1, "b--", zorder=4, transform=geo)

    for i in range(len(names)):
        # Create text box with semi-opaque background
        bbox = dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="none")

        # Add label with text box and transform for projection
        ax.text(
            x[i],
            y[i],
            names[i],
            va="top",
            family="monospace",
            weight="bold",
            bbox=bbox,
            transform=geo,
        )

    # Add scale bar
    scale_bar(ax, (0.1, 0.1), 50)

    if set_global:
        ax.set_global()
    if setup_name != "":
        setup_name = f"{setup_name}_"

    if not os.path.exists("plots"):
        os.makedirs("plots")
    fn = f"plots/{setup_name}station_map.pdf"
    plt.savefig(fn, bbox_inches="tight")
    print(f"done writing {fn}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate teleseismic station map")
    parser.add_argument("config_file", help="config file describing event and stations")
    parser.add_argument(
        "--plot_all_station_file",
        action="store_true",
    )

    args = parser.parse_args()

    config = configparser.ConfigParser()
    assert os.path.isfile(args.config_file), f"{args.config_file} not found"
    config.read(args.config_file)

    setup_name = config.get("GENERAL", "setup_name")
    lSourceFiles = config.get("GENERAL", "source_files").split(",")
    client_name = config.get("GENERAL", "client")
    onset = config.get("GENERAL", "onset")
    event_lon = config.getfloat("GENERAL", "hypo_lon")
    event_lat = config.getfloat("GENERAL", "hypo_lat")
    event = {"longitude": event_lon, "latitude": event_lat, "onset": onset}
    station_codes = config.get("GENERAL", "stations").split(",")
    station_codes = [x.strip() for x in station_codes]

    software = config.get("GENERAL", "software").split(",")
    set_global = "axitra" not in software
    station_file = config.get("GENERAL", "station_file", fallback=None)
    path_observations = config.get("GENERAL", "path_observations")

    if args.plot_all_station_file:
        df = pd.read_csv(station_file)
        df.rename(columns={"lon": "longitude", "lat": "latitude"}, inplace=True)
    else:
        df = retrieve_coordinates(client_name, event, station_codes)
    print(df)

    faultfname = config.get("GENERAL", "fault_filename", fallback=None)

    if faultfname:
        polygons = compute_shapely_polygon(faultfname)
        projection = config.get("GENERAL", "projection")
        fault_info = {}
        fault_info["projection"] = projection
        fault_info["polygons"] = polygons
    else:
        fault_info = None

    generate_station_map(df, event, set_global, setup_name, fault_info=fault_info)
