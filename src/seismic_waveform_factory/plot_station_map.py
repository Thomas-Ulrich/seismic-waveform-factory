#!/usr/bin/env python3
import argparse
import glob
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from obspy import UTCDateTime, read_inventory
from obspy.core.inventory import Inventory
from pyproj import Transformer
from geographiclib.geodesic import Geodesic

from seismic_waveform_factory.config.loader import ConfigLoader
from seismic_waveform_factory.config.schema import CONFIG_SCHEMA
from seismic_waveform_factory.config.utils import categorize_stations_by_scale
from seismic_waveform_factory.fault.fault_processing import compute_shapely_polygon
from seismic_waveform_factory.geo.utils import add_distance_backazimuth_to_df
from seismic_waveform_factory.utils.scalebar import scale_bar
from seismic_waveform_factory.waveform.retrieve import initialize_client


def add_great_circle_radius(ax, event, radius_deg=30):
    geod = Geodesic.WGS84
    lons_circle = []
    lats_circle = []
    for az in np.linspace(0, 360, 361):
        g = geod.ArcDirect(event["lat"], event["lon"], az, radius_deg)
        lons_circle.append(g["lon2"])
        lats_circle.append(g["lat2"])
    ax.plot(
        lons_circle,
        lats_circle,
        transform=ccrs.Geodetic(),
        color="grey",
        linestyle="--",
        linewidth=0.5,
        zorder=2,
    )
    # --- Annotate the radius at the bottom of the circle (azimuth 180°) ---
    g_bottom = geod.ArcDirect(event["lat"], event["lon"], 180, radius_deg)
    lon_text, lat_text = g_bottom["lon2"], g_bottom["lat2"]

    ax.text(
        lon_text,
        lat_text - 0.5,  # small offset southward for readability
        f"{radius_deg}°",
        color="k",
        ha="center",
        va="top",
        transform=ccrs.Geodetic(),
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1.5),
    )

    return ax


def retrieve_coordinates(cfg, station_codes):
    path_observations = cfg["general"]["path_observations"]
    os.makedirs(path_observations, exist_ok=True)

    rows = []
    fn_inventories = glob.glob(os.path.join(path_observations, "*.xml"))
    for fn_inventory in fn_inventories:
        inventory = read_inventory(fn_inventory)
        for network in inventory:
            for station in network:
                code = f"{network.code}.{station.code}"
                if code in station_codes:
                    rows.append(
                        {
                            "network": network.code,
                            "station": station.code,
                            "longitude": station.longitude,
                            "latitude": station.latitude,
                        }
                    )

    df = pd.DataFrame(rows, columns=["network", "station", "longitude", "latitude"])
    df["code"] = df["network"] + "." + df["station"]
    print(df)

    missing = [code for code in station_codes if code not in df["code"].values]
    if missing:
        print(f"{missing} locations not found in inventory files:\n{fn_inventories}")
        client_name = cfg["general"]["client"]
        client = initialize_client(client_name)
        event_time = UTCDateTime(cfg["general"]["hypocenter"]["onset"])
        networks = set([entry.split(".")[0] for entry in missing])
        inv = Inventory()
        for net in networks:
            print(f"retrieving station data from network {net}")
            inventory = client.get_stations(
                network=net,
                level="station",
                starttime=event_time,
                endtime=event_time + 100,
                includeavailability=True,
            )
            fn_inventory = f"{path_observations}/inv_{client_name}.xml"
            inventory.write(fn_inventory, format="STATIONXML")
            inv.extend(inventory)

        network_station_pairs = [code.split(".") for code in missing]
        extracted_data = [
            (network.code, station.code, station.longitude, station.latitude)
            for network in inv.networks
            for station in network.stations
            if [network.code, station.code] in network_station_pairs
        ]
        print(extracted_data)
        rows = []
        for vals in extracted_data:
            net, sta, lon, lat = vals
            rows.append(
                {
                    "network": net,
                    "station": sta,
                    "longitude": lon,
                    "latitude": lat,
                }
            )
        dfm = pd.DataFrame(
            rows, columns=["network", "station", "longitude", "latitude"]
        )
        dfm["code"] = dfm["network"] + "." + dfm["station"]
        df = pd.concat([df, dfm], ignore_index=True)

    df = add_distance_backazimuth_to_df(df, cfg["general"]["hypocenter"])
    return df


def generate_station_map(df, cfg, set_global=False):
    plt.figure(figsize=(6, 6))
    event = cfg["general"]["hypocenter"]

    faultfname = cfg["general"]["fault_file"]
    if faultfname:
        polygons = compute_shapely_polygon(faultfname)
        projection = cfg["general"]["projection"]
        assert projection, "projection is need to plot the fault"
        fault_info = {}
        fault_info["projection"] = projection
        fault_info["polygons"] = polygons
    else:
        fault_info = None

    if set_global:
        projection = ccrs.Orthographic(
            central_longitude=event["lon"], central_latitude=event["lat"]
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
        x1, x2 = min(event["lon"], x.min()), max(event["lon"], x.max())
        y1, y2 = min(event["lat"], y.min()), max(event["lat"], y.max())
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
        event["lon"],
        event["lat"],
        200,
        color="orange",
        marker="*",
        edgecolor="k",
        zorder=3,
        transform=geo,
    )

    faultfname = cfg["general"]["fault_file"]
    if faultfname:
        polygons = compute_shapely_polygon(faultfname)
        projection = cfg["general"]["projection"]
        assert projection, "projection is need to plot the fault"

        transformer = Transformer.from_crs(projection, "epsg:4326", always_xy=True)
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
        map_type = "teleseismic"
        ax = add_great_circle_radius(ax, event, radius_deg=30)
        ax = add_great_circle_radius(ax, event, radius_deg=60)
    else:
        map_type = "regional"
        # Add scale bar
        scale_bar(ax, (0.1, 0.1), 50)

    setup_name = cfg["general"]["setup_name"]
    if setup_name:
        setup_name += "_"
    if not os.path.exists("plots"):
        os.makedirs("plots")
    fn = f"plots/{setup_name}station_map_{map_type}.pdf"
    plt.savefig(fn, bbox_inches="tight")
    print(f"done writing {fn}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate station map")
    parser.add_argument("config_file", help="config file describing event and stations")
    parser.add_argument(
        "--plot_all_station_file",
        action="store_true",
    )
    args = parser.parse_args()
    cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)

    station_codes = categorize_stations_by_scale(cfg)

    if args.plot_all_station_file:
        station_file = cfg["general"]["station_file"]
        assert station_file, "station_file is needed with plot_all_station_file"
        df = pd.read_csv(station_file)
        df.rename(columns={"lon": "longitude", "lat": "latitude"}, inplace=True)
        set_global = False
        generate_station_map(df, cfg, set_global)
    else:
        for v in station_codes.keys():
            set_global = True if v == "global" else False
            stations = station_codes[v]
            print(v, stations)
            if stations:
                df = retrieve_coordinates(cfg, stations)
                generate_station_map(df, cfg, set_global)
