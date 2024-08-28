from obspy import read
from lxml.etree import XMLSyntaxError
from obspy.clients.fdsn.header import FDSNNoDataException, FDSNException
from obspy.clients.fdsn import Client, RoutingClient
from requests.exceptions import ConnectionError
from obspy.core.inventory import Inventory
from obspy.core.util.obspy_types import ObsPyException
from obspy import read_inventory
import gzip
import os


def load_cached_station_data(network, stations, level, cache_dir):
    inv = Inventory()
    stations_not_cached = []
    for station in stations:
        cache_file = os.path.join(cache_dir, f"{network}_{station}_{level}.xml")
        if os.path.exists(cache_file):
            print(f"Loading cached data for {network}.{station}")
            station_inv = read_inventory(cache_file)
            inv.extend(station_inv)
        else:
            print(f"No cached data found for {network}.{station}")
            stations_not_cached.append(station)
    return inv, stations_not_cached


def save_station_data(inventory, level, cache_dir):
    # Save inventory station-wise regardless of retrieval method
    for net in inventory:
        for sta in net:
            station_inv = Inventory(networks=[net.copy()])
            station_inv[0].stations = [sta]
            station_cache_file = os.path.join(
                cache_dir, f"{net.code}_{sta.code}_{level}.xml"
            )
            station_inv.write(station_cache_file, format="STATIONXML")
            print(f"Saved data for station {net.code}.{sta.code}")


def get_station_data(client, network, stations, level, t1, network_wise=True):
    exceptions_to_catch = (
        FDSNException,
        FDSNNoDataException,
        XMLSyntaxError,
        gzip.BadGzipFile,
        ConnectionError,
    )
    cache_dir = "observations"
    os.makedirs(cache_dir, exist_ok=True)
    inv, stations_not_cached = load_cached_station_data(
        network, stations, level, cache_dir
    )
    if not stations_not_cached:
        return inv
    if network_wise:
        station_param = [",".join(stations_not_cached)]
    else:
        station_param = stations_not_cached

    for station in station_param:
        retry_message = f"network {network} at station(s) {station}"
        max_retries = 5
        for retry_count in range(max_retries):
            try:
                if retry_count == 0:
                    print(f"Getting {level}s for {retry_message}...")
                inventory = client.get_stations(
                    network=network,
                    station=station,
                    level=level,
                    starttime=t1,
                )

                inv.extend(inventory)
                save_station_data(inventory, level, cache_dir)
                break
            except exceptions_to_catch as e:
                if retry_count == max_retries - 1:
                    print(f"Max retry count reached for {retry_message}. Skipping.")
                    print(f"Failed getting {level} for {retry_message}")
                    print(f"Error: {e.__class__.__name__}")
                else:
                    print(
                        f"Error occurred in get_station for {retry_message}: {e.__class__.__name__}"
                    )
                continue
    return inv


def get_waveforms(
    client, network, station, selected_band, t1, t2, is_routing_client=False
):
    max_retries = 5
    exceptions_to_catch = (
        FDSNException,
        FDSNNoDataException,
        XMLSyntaxError,
        gzip.BadGzipFile,
        ConnectionError,
    )
    st_obs0 = False

    if is_routing_client:
        kwargs = {}
    else:
        kwargs = {"attach_response": True}

    for retry_count in range(max_retries):
        try:
            st_obs0 = client.get_waveforms(
                network=network,
                station=station.code,
                location="*",
                channel=f"{selected_band}*",
                starttime=t1,
                endtime=t2,
                **kwargs,
            )
            if not st_obs0:
                print(f"Got empty stream for {network} {station.code}")
            break
        except exceptions_to_catch as e:
            if retry_count == max_retries - 1:
                print(f"Max retry count reached for {station.code}. Skipping.")
                print(f"Error: {e.__class__.__name__}")
            else:
                print(
                    f"Error occurred in get_waveforms for {station.code}: {e.__class__.__name__}"
                )
            continue

    if not st_obs0:
        print(f"No waveform available for {station.code}")

    return st_obs0


def get_pre_filt(selected_band):
    if selected_band == "L":
        return [0.00033, 0.001, 0.1, 0.3]
    else:
        return [0.001, 0.005, 45, 50]


def select_band(channels):
    priorities = ["H", "B", "E", "M", "L"]
    for band in priorities:
        if any(channel.startswith(band) for channel in channels):
            if any(channel.startswith(band + "H") for channel in channels):
                return band + "H"
            elif any(channel.startswith(band + "N") for channel in channels):
                return band + "N"
            else:
                return band
    return None


def retrieve_waveforms(
    network_station, client_name, kind_vd, path_observations, t1, t2
):
    if client_name in ["eida-routing", "iris-federator"]:
        client = RoutingClient(client_name)
        is_routing_client = True
        level = "response"
    else:
        client = Client(client_name)
        is_routing_client = False
        level = "channel"

    os.makedirs(path_observations, exist_ok=True)
    retrieved_waveforms = {}
    for network, stations in network_station.items():
        for station in stations.copy():
            fname = f"{network}.{station}_{kind_vd}_{t1.date}.mseed"
            fullfname = os.path.join(path_observations, fname)
            if os.path.exists(fullfname):
                print(f"reading the data from {fullfname}")
                code = f"{network}.{station}"
                retrieved_waveforms[code] = read(fullfname)
                stations.remove(station)
        if not stations:
            continue

        if level == "channel":
            inventory = get_station_data(
                client, network, stations, level, t1, network_wise=True
            )
        else:
            inventory = get_station_data(
                client, network, stations, level, t1, network_wise=False
            )

        if len(inventory) == 0:
            print(f"could not get {level} for {network}")
            continue

        for station in [sta for net in inventory for sta in net]:
            code = f"{network}.{station.code}"
            fname = f"{code}_{kind_vd}_{t1.date}.mseed"
            fullfname = os.path.join(path_observations, fname)

            print(f"requesting data for {network}.{station.code}")
            channels = [channel.code for channel in station]

            def has_band(x):
                return any(channel.startswith(x) for channel in channels)

            selected_band = select_band(channels)
            if not selected_band:
                print(f"{station.code} does not have the expected channels: {channels}")
                continue
            st_obs0 = get_waveforms(
                client, network, station, selected_band, t1, t2, is_routing_client
            )

            if not st_obs0:
                continue

            st_obs0.rotate(method="->ZNE", inventory=inventory)

            # define a filter band to prevent amplifying noise during the deconvolution
            # pre_filt = [0.00033, 0.001, 1.0, 3.0]
            output_dic = {
                "acceleration": "ACC",
                "velocity": "VEL",
                "displacement": "DISP",
            }

            pre_filt = get_pre_filt(selected_band)
            exceptions_to_catch = (ValueError, ObsPyException)
            try:
                st_obs0.remove_response(
                    output=output_dic[kind_vd],
                    pre_filt=pre_filt,
                    # todo: use water_level if kind_vd == instrument measured quantity (see warning in)
                    # https://docs.obspy.org/master/packages/autogen/obspy.core.trace.Trace.remove_response.html
                    water_level=None,
                    zero_mean=True,
                    taper=True,
                    taper_fraction=0.05,
                    inventory=inventory,
                )
            except exceptions_to_catch as e:
                # In theory this should not happen, but it does...
                # ValueError: No response information found. Use `inventory` parameter to specify an inventory with response information.
                # or
                # obspy.core.util.obspy_types.ObsPyException: Can not use evalresp on response with no response stages.
                print(
                    f"Error in st_obs0.remove_response  at station {code}: {e.__class__.__name__}"
                )
                continue
            st_obs0.write(fullfname, format="MSEED")
            retrieved_waveforms[code] = st_obs0
            print(f"done writing {fullfname}")
    return retrieved_waveforms
