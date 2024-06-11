import os
import numpy as np
from obspy import read, Trace
from obspy.clients.fdsn.header import FDSNException
import xml.etree.ElementTree as ET
import pyproj
import matplotlib.pyplot as plt
import glob
import pandas as pd
from pyproj import Transformer
from obspy import Stream

import glob
import numpy as np

def read_seissol_receiver_file(output_path, idst, coords_only=False):
    """
    Read SeisSol receiver seismogram data.

    Args:
        output_path (str): The path to the output directory, including the prefix.
        idst (int): The receiver ID.
        coords_only (bool, optional): If True, return only the receiver coordinates. Defaults to False.

    Returns:
        tuple or list: If `coords_only` is False, returns a tuple containing the receiver coordinates,
                       variable list, and seismogram data.
                       If `coords_only` is True, returns a list containing the receiver coordinates.
    """
    # Construct the file pattern for the receiver data
    file_pattern = f"{output_path}-receiver-{idst:05d}*"
    files = glob.glob(file_pattern)

    if not files:
        print(f"No files found matching the pattern: {file_pattern}")
        return None

    if len(files) > 1:
        print(f"Warning: Multiple files match the pattern {file_pattern}. Using the first one: {files[0]}")

    with open(files[0], 'r') as file:
        # Skip the first line
        file.readline()

        # Read and process the variable list
        variable_list = file.readline()[11:].split(",")
        variable_list = np.array([var.strip().strip('"') for var in variable_list])

        # Read the receiver coordinates
        xsyn = float(file.readline().split()[2])
        ysyn = float(file.readline().split()[2])
        zsyn = float(file.readline().split()[2])

        if coords_only:
            return [xsyn, ysyn, zsyn]

        # Read the seismogram data
        seismogram = np.loadtxt(file)

    return ([xsyn, ysyn, zsyn], variable_list, seismogram)

def get_station_code_from_coordinates(inventory, lonlatdepth, eps=5e-3):
    """
    Find the station code (e.g., KO.FOCM, HL.TNSA) from the given coordinates.

    Args:
        inventory (obspy.core.inventory.Inventory): The inventory containing station information.
        lonlatdepth (tuple or list): A tuple or list containing longitude, latitude, and depth values.
        eps (float, optional): The maximum allowable difference between coordinates (default: 5e-3).

    Returns:
        str or None: The station code if found, or None if not found.
    """
    target_lon, target_lat, _ = lonlatdepth

    for network in inventory:
        for station in network:
            station_lon = station.longitude
            station_lat = station.latitude

            if abs(station_lon - target_lon) < eps and abs(station_lat - target_lat) < eps:
                station_code = f"{network.code}.{station.code}"
                print(f"Station code: {station_code}, Longitude: {station_lon}, Latitude: {station_lat}")
                return station_code
    return None


def stream_from_seissol_data(network_code, station_code, variable_list, synth, starttime):
    """
    Load SeisSol receiver data into an ObsPy Stream object.

    Args:
        network_code (str): The network code.
        station_code (str): The station code.
        variable_list (list or np.ndarray): A list or array of variable names.
        synth (np.ndarray): The seismogram data.
        starttime (obspy.UTCDateTime): The start time of the seismogram.

    Returns:
        obspy.Stream: A Stream object containing the seismogram data.
    """
    st_syn = Stream()
    xyz = "ENZ"
    uvw = ["u", "v", "w"] if "u" in variable_list else [f"v{i}" for i in range(1, 4)]

    for i, channel in enumerate(xyz):
        try:
            j = variable_list.tolist().index(uvw[i])
        except ValueError:
            print(f"Variable '{uvw[i]}' not found in the variable list: {variable_list}")
            continue

        tr = Trace()
        tr.stats.station = station_code
        tr.stats.network = network_code
        tr.stats.channel = channel
        tr.data = synth[:, j]
        tr.stats.delta = synth[1, 0] - synth[0, 0]
        tr.stats.starttime = starttime
        st_syn.append(tr)

    return st_syn


def create_zero_stream(network_code, station_code, starttime, delta=1.0, npts=2):
    """
    Create an ObsPy Stream object with zero-valued trace data.

    Args:
        network_code (str): The network code.
        station_code (str): The station code.
        starttime (obspy.UTCDateTime): The start time of the trace.
        delta (float, optional): The time delta (sample spacing) in seconds. Default is 1.0.
        npts (int, optional): The number of data points for each trace. Default is 2.

    Returns:
        obspy.Stream: A Stream object containing three zero-valued traces (E, N, Z components).
    """
    st_syn = Stream()
    xyz = "ENZ"

    for channel in xyz:
        tr = Trace()
        tr.stats.station = station_code
        tr.stats.network = network_code
        tr.stats.channel = channel
        tr.data = np.zeros(npts)
        tr.stats.delta = delta
        tr.stats.starttime = starttime
        st_syn.append(tr)

    return st_syn

def compile_inv_lut_gm(folder_prefix, projection, inventory):
    """
    Compile an inverse lookup table for station codes and receiver IDs.

    Args:
        folder_prefix (str): The prefix of the folder containing the receiver data.
        projection (str or pyproj.Proj): The projection to use for coordinate transformations.
        inventory (obspy.core.inventory.Inventory or None, optional): The inventory containing station information.

    Returns:
        dict: A dictionary mapping station codes to receiver IDs for the specified stations.
    """
    transformer = Transformer.from_crs(projection, "epsg:4326", always_xy=True)

    print("Compiling lookup table...")
    station_lookup_table = {}
    id_not_found = []

    file_pattern = f"{folder_prefix}-receiver-*"
    receiver_files = glob.glob(file_pattern)

    for fn in receiver_files:
        id_station = int(fn.split(f"{folder_prefix}-receiver-")[1].split("-")[0])

        # Load SeisSol receiver coordinates
        xyzs = read_seissol_receiver_file(folder_prefix, id_station, coords_only=True)
        lonlatdepth = transformer.transform(xyzs[0], xyzs[1], xyzs[2])

        station_code = get_station_code_from_coordinates(inventory, lonlatdepth)
        if station_code:
            station_lookup_table[id_station] = station_code
        else:
            id_not_found.append(id_station)

    if id_not_found:
        print(f"Warning: No station codes found for receiver IDs: {id_not_found}")

    inv_station_lookup_table = {v: k for k, v in station_lookup_table.items()}
    return inv_station_lookup_table

def collect_seissol_synthetics(seissol_outputs, list_inventory, projection, t1):
    """
    Collect synthetic seismograms from SeisSol outputs.

    Args:
        seissol_outputs (list): A list of paths to SeisSol output directories.
        list_inventory (list): A list of ObsPy inventory objects.
        projection (str or pyproj.Proj): The projection to use for coordinate transformations.
        t1 (obspy.UTCDateTime): The start time of the seismograms.

    Returns:
        list: A list of ObsPy Stream objects containing synthetic seismograms.
    """
    combined_inventory = list_inventory[0]
    for inv in list_inventory[1:]:
        combined_inventory += inv

    inv_station_lookup_tables = []
    for seissol_output_path in seissol_outputs:
        inv_station_lookup_tables.append(
            compile_inv_lut_gm(seissol_output_path, projection, combined_inventory)
        )

    list_synthetics = []
    for idsyn, seissol_output in enumerate(seissol_outputs):
        syn_st = Stream()
        for net in combined_inventory:
            for station in net:
                station_code = f"{net.code}.{station.code}"
                if station_code in inv_station_lookup_tables[idsyn]:
                    id_station = inv_station_lookup_tables[idsyn][station_code]
                    xyzs, variablelist, synth = read_seissol_receiver_file(
                        seissol_output, id_station
                    )
                    syn_st += stream_from_seissol_data(
                        net.code, station.code, variablelist, synth, t1
                    )
                else:
                    print("Station {station_code} not found in SeisSol receivers")
                    syn_st += create_zero_stream(net.code, station.code, t1)
        list_synthetics.append(syn_st)
    return list_synthetics



