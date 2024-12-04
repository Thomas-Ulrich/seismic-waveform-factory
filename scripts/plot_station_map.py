#!/usr/bin/env python3
from obspy import UTCDateTime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse
import configparser
from obspy.clients.fdsn import Client, RoutingClient
import os
from fault_processing import compute_shapely_polygon
from retrieve_waveforms import get_station_data
from waveform_figure_utils import initialize_client


def retrieve_coordinates(client_name, event, station_codes):
    client = initialize_client(client_name)

    event_time = UTCDateTime(event["onset"])

    df = pd.DataFrame(columns=["network", "station", "longitude", "latitude"])

    for ins, netStaCode in enumerate(station_codes):
        listNetStaCode = netStaCode.split(".")
        if len(listNetStaCode) == 1:
            station = netStaCode
            network = "*"
        else:
            network, station = listNetStaCode

        inventory = get_station_data(
            client, network, [station], "station", event_time, network_wise=False
        )
        new_row = {
            "network": inventory[0].code,
            "station": inventory[0][0].code,
            "longitude": inventory[0][0].longitude,
            "latitude": inventory[0][0].latitude,
        }
        df.loc[len(df)] = new_row
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
        ax.set_extent([x1 - dx, x2 + dx, y1 - dy, y2 + dy], ccrs.PlateCarree())
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
    if set_global:
        ax.set_global()
    if setup_name != "":
        setup_name = f"{setup_name}_"

    if not os.path.exists("plots"):
        os.makedirs("plots")
    fn = f"plots/{setup_name}station_map.svg"
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
    software = config.get("GENERAL", "software").split(",")
    set_global = "axitra" not in software
    station_file = config.get("GENERAL", "station_file", fallback=None)

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
