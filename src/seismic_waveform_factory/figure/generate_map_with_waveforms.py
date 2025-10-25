#!/usr/bin/env python3
import glob
import os

# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, read_inventory
from obspy.core import UTCDateTime

from seismic_waveform_factory.config.loader import ConfigLoader
from seismic_waveform_factory.config.schema import CONFIG_SCHEMA


def main(args):
    cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)

    # Extract station information
    panel = args.plot_panel
    wf_plot = cfg["waveform_plots"][panel]
    station_list = wf_plot["stations"]
    kind_vd = wf_plot["kind"]
    print(station_list)

    event = cfg["general"]["hypocenter"]
    hypo_lat = event["lat"]
    hypo_lon = event["lon"]
    origin_time = UTCDateTime(event["onset"])

    filter_fmin = 1.0 / wf_plot["filter_tmax"]
    filter_fmax = 1.0 / wf_plot["filter_tmin"]

    # Set paths
    waveform_folder = cfg["general"]["path_observations"]

    # Initialize station dictionary
    stations = {}

    # Find waveform and response files
    for station in station_list:
        net, sta = station.split(".")

        # Find waveform file
        pattern = f"{waveform_folder}/{net}.{sta}_*H_{kind_vd}_*.mseed"
        waveform_file = glob.glob(pattern)

        if not waveform_file:
            print(f"Warning: No waveform found for {station} with pattern {pattern}")
            continue
        waveform_file = waveform_file[0]

        # Find response file
        response_file = glob.glob(f"{waveform_folder}/{net}_{sta}_response.xml")
        if not response_file:
            print(f"Warning: No response found for {station}")
            continue
        response_file = response_file[0]

        # Store in dictionary
        stations[station] = {"waveform": waveform_file, "response": response_file}

    lats = [hypo_lat]
    lons = [hypo_lon]
    # Plot stations with waveforms
    for station, info in stations.items():
        net, sta = station.split(".")
        inv = read_inventory(info["response"])
        y, x = (
            inv.networks[0].stations[0].latitude,
            inv.networks[0].stations[0].longitude,
        )
        lats.append(y)
        lons.append(x)

    # Calculate map boundaries with padding
    lat_min = min(lats) - 0.5
    lat_max = max(lats) + 0.5
    lon_min = min(lons) - 0.5
    lon_max = max(lons) + 0.8

    # Create figure and axes
    fig, ax = plt.subplots(
        figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()}
    )

    # Set map boundaries
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add geographic features
    ax.add_feature(cfeature.LAND, facecolor="lightgray")
    plt.plot()
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    # ax.add_feature(cfeature.STATES, linestyle=':')

    # Add gridlines
    gl = ax.gridlines(
        draw_labels=True,
        # xlocs=np.arange(round_min_down(lon_min), round_max_up(lon_max), 0.5),
        xlocs=np.arange(-180, 180, 0.5),
        ylocs=np.arange(-90, 90, 0.5),
        linewidth=0.5,
        color="gray",
        alpha=0.5,
        linestyle="--",
        dms=False,
        x_inline=False,
        y_inline=False,
    )
    print(np.arange(lon_min, lon_max, 0.5))
    gl.top_labels = False
    gl.right_labels = False

    ax.scatter(
        hypo_lon,
        hypo_lat,
        marker="*",
        color="yellow",
        edgecolor="black",
        s=200,
        zorder=3,
    )
    plt.plot()

    # Plot stations with waveforms
    for station, info in stations.items():
        net, sta = station.split(".")

        # Read waveform
        st = read(info["waveform"])

        # Read response
        inv = read_inventory(info["response"])

        # Remove instrument response
        st.remove_response(inventory=inv, output="ACC")  # Output in acceleration

        st.filter(
            "bandpass",
            freqmin=filter_fmin,
            freqmax=filter_fmax,
            corners=4,
            zerophase=True,
        )
        st.trim(starttime=origin_time, endtime=origin_time + 100)

        # Normalize waveform for plotting
        tr = st[0]  # First trace
        waveform = tr.data / max(abs(tr.data)) * 0.3  # Scale factor
        # Compute time relative to origin_time using reftime
        time = tr.times(reftime=origin_time)

        # Get station coordinates from response file
        y, x = (
            inv.networks[0].stations[0].latitude,
            inv.networks[0].stations[0].longitude,
        )

        # Plot station location
        ax.scatter(x, y, marker="^", color="red", edgecolor="black", s=100, zorder=3)
        ax.text(x, y + 0.1, station, ha="center", color="darkblue", zorder=4)

        # Plot waveform near station
        ax.plot(x + time * 0.005, y + waveform * 0.2, "k", lw=1)
    fname = "map_with_waveform.png"
    plt.savefig(fname)
    print(f"done generating {fname}")
    full_path = os.path.abspath(fname)
    print(f"full path: {full_path}")

    # Show the plot
    # plt.show()
