import pandas as pd
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth
from obspy.clients.fdsn import Client, RoutingClient
import functools as ft
from lxml.etree import XMLSyntaxError
import gzip
import pickle


def compile_list_inventories(client_name, station_codes, t1):
    max_retries = 5
    retry_count = 0
    while retry_count < max_retries:
        try:
            if client_name == "eida-routing":
                c = RoutingClient(client_name)
            else:
                c = Client(client_name)
            break
        except UnicodeDecodeError as e:
            print("excepting UnicodeDecodeError in Client", e)
            retry_count += 1
            if retry_count == max_retries:
                raise UnicodeDecodeError("Max retry count reached")

    list_inventory = []
    for ins, netStaCode in enumerate(station_codes):
        listNetStaCode = netStaCode.split(".")
        if len(listNetStaCode) == 1:
            station = netStaCode
            network = "*"
        else:
            network, station = listNetStaCode
        max_retries = 5
        retry_count = 0
        while retry_count < max_retries:
            try:
                if retry_count == 0:
                    print(f"getting channels for network {network}...")
                inventory = c.get_stations(
                    network=network, station=station, level="channel", starttime=t1
                )
                break
            except XMLSyntaxError as e:
                print(
                    f"excepting XML Syntax Error in get_stations for network {network}",
                    e,
                )
                retry_count += 1
                if retry_count == max_retries:
                    raise XMLSyntaxError("Max retry count reached")
            except gzip.BadGzipFile as e:
                print(
                    f"excepting gzip.BadGzipFile in get_stations for network {network}",
                    e,
                )
                retry_count += 1
                if retry_count == max_retries:
                    raise gzip.BadGzipFile("Max retry count reached")
            except ValueError as e:
                # ValueError: The current client does not have a station service.
                print(f"excepting ValueError in get_stations for network {network}", e)
                retry_count += 1
                if retry_count == max_retries:
                    raise ValueError("Max retry count reached")
        # print(inventory)
        list_inventory.append(inventory)

    return list_inventory


def reorder_station_using_azimuth(list_inventory, hypo_lon, hypo_lat):
    # Reorder station based on azimuth
    d_azimuth = {}
    for ins, inv in enumerate(list_inventory):
        sta = inv[0][0]
        d_azimuth[ins] = gps2dist_azimuth(
            lat1=sta.latitude, lon1=sta.longitude, lat2=hypo_lat, lon2=hypo_lon
        )[2]
    ordered_azimuth = {
        k: v for k, v in sorted(d_azimuth.items(), key=lambda item: item[1])
    }
    list_inventory2 = list_inventory.copy()
    list_inventory = []
    for key in ordered_azimuth:
        list_inventory.append(list_inventory2[key])
    return list_inventory


def estimate_travel_time(source_depth_in_km, distance_in_degree, station, phase="P"):
    taupModel = "ak135"
    model = TauPyModel(model=taupModel)
    tP = model.get_travel_times(
        source_depth_in_km=source_depth_in_km,
        distance_in_degree=distance_in_degree,
        phase_list=[phase],
    )
    if not tP:
        print(f"no P wave at station {station}")
        tP = 0.0
    else:
        tP = tP[0].time
    return tP


def merge_gof_dfs(Pwave, SHwave, surface_waves):
    gofall_dfs = []
    if Pwave.enabled:
        gofall_dfs.append(Pwave.gof_df)

    if SHwave.enabled:
        gofall_dfs.append(SHwave.gof_df)

    if surface_waves.enabled:
        gofall_dfs.append(surface_waves.gof_df)

    df_final = ft.reduce(
        lambda left, right: pd.merge(
            left, right, on="station", suffixes=("", "_remove")
        ),
        gofall_dfs,
    )
    df_final.drop([i for i in df_final.columns if "remove" in i], axis=1, inplace=True)
    return df_final
