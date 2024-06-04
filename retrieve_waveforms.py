from obspy import read
from lxml.etree import XMLSyntaxError
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.clients.fdsn import Client, RoutingClient
from obspy.core.inventory import Inventory
import gzip
import os


def get_level_station_wise(client, network, stations, level, t1):
    inv = Inventory()
    exceptions_to_catch = (FDSNNoDataException, XMLSyntaxError, gzip.BadGzipFile)

    for station in stations:
        max_retries = 5
        for retry_count in range(max_retries):
            try:
                if retry_count == 0:
                    print(
                        f"Getting {level}s for network {network} and station {station}..."
                    )
                inventory = client.get_stations(
                    network=network,
                    station=station,
                    level=level,
                    starttime=t1,
                )
                inv.extend(inventory)
                break
            except exceptions_to_catch as e:
                if retry_count == max_retries - 1:
                    print(f"Max retry count reached for {network}. Skipping.")
                    print(f"Failed getting {level} for {station}")
                    print(f"Error: {e.__class__.__name__}")
                else:
                    print(
                        f"Error occurred in get_station for network {network} at station {station}: {e.__class__.__name__}"
                    )
                continue

    return inv


def get_level_network_wise(client, network, stations, level, t1):
    max_retries = 5
    inv = Inventory()
    exceptions_to_catch = (FDSNNoDataException, XMLSyntaxError, gzip.BadGzipFile)

    for retry_count in range(max_retries):
        try:
            if retry_count == 0:
                print(f"Getting {level}s for network {network}...")
            inv = client.get_stations(
                network=network,
                station=",".join(stations),
                level=level,
                starttime=t1,
            )
            break
        except exceptions_to_catch as e:
            if retry_count == max_retries - 1:
                print(f"Max retry count reached for {network}. Skipping.")
                print(f"Error: {e.__class__.__name__}")
            else:
                print(
                    f"Error occurred in get_station for network {network}: {e.__class__.__name__}"
                )
            continue

    return inv


def get_waveforms(
    client, network, station, selected_band, t1, t2, is_routing_client=False
):
    max_retries = 5
    exceptions_to_catch = (FDSNNoDataException, XMLSyntaxError)
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


def retrieve_waveforms(
    network_station, client_name, kind_vd, path_observations, t1, t2
):
    if client_name == "eida-routing":
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

        if "level" == "channel":
            inventory = get_level_network_wise(client, network, stations, level, t1)
        else:
            inventory = get_level_station_wise(client, network, stations, level, t1)

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

            if has_band("H"):
                selected_band = "H"
            elif has_band("B"):
                selected_band = "B"
            elif has_band("M"):
                selected_band = "M"
            elif has_band("L"):
                selected_band = "L"
            else:
                print(f"{station.code} has not the expected channels: {channels}")
                continue

            # now considering the second letter
            if has_band(selected_band + "H"):
                selected_band += "H"
            elif has_band(selected_band + "N"):
                selected_band += "N"
            else:
                print(f"{station.code} has not the expected channels: {channels}")
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
            pre_filt = [0.001, 0.005, 45, 50]
            if selected_band == "L":
                pre_filt = [0.00033, 0.001, 0.1, 0.3]
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
            except ValueError:
                # In theory this should not happen, but it does...
                # ValueError: No response information found. Use `inventory` parameter to specify an inventory with response information.
                print("ValueError: when removing_response at station {code}")
                continue
            st_obs0.write(fullfname, format="MSEED")
            retrieved_waveforms[code] = st_obs0
            print(f"done writing {fullfname}")
    return retrieved_waveforms
