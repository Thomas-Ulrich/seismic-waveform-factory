import numpy as np
import configparser
import glob
import matplotlib.pyplot as plt
from obspy import read, read_inventory
from obspy.core import UTCDateTime

# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Load the configuration file
config = configparser.ConfigParser()
config.read("waveforms_config.ini")

# Extract station information
station_list = [st.strip() for st in config["GENERAL"]["stations"].split(",")]

hypo_lat = float(config["GENERAL"]["hypo_lat"])
hypo_lon = float(config["GENERAL"]["hypo_lon"])
hypo_depth = float(config["GENERAL"]["hypo_depth_in_km"])
origin_time = UTCDateTime(config["GENERAL"]["onset"])

filter_fmin = 1.0 / config.getfloat("P_WAVE", "filter_tmax")
filter_fmax = 1.0 / config.getfloat("P_WAVE", "filter_tmin")

# Set paths
waveform_folder = "observations/"

# Initialize station dictionary
stations = {}

# Find waveform and response files
for station in station_list:
    net, sta = station.split(".")

    # Find waveform file
    waveform_file = glob.glob(f"{waveform_folder}{net}.{sta}_HN_acceleration_*.mseed")
    if not waveform_file:
        print(f"Warning: No waveform found for {station}")
        continue
    waveform_file = waveform_file[0]

    # Find response file
    response_file = glob.glob(f"{waveform_folder}{net}_{sta}_response.xml")
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
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()})

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
plt.savefig("map.png")
# Show the plot
# plt.show()
