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
from concurrent.futures import ThreadPoolExecutor, as_completed


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
        level = "response"

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
            try:
                st_obs0.rotate(method="->ZNE", inventory=inventory)
            except ValueError as e:
                # get rid of this rare error:
                # raise ValueError("The given directions are not linearly independent, "
                # ValueError: The given directions are not linearly independent,
                # at least within numerical precision. Determinant of the base change matrix: 0
                print(
                    f"Error in st_obs0.rotate  at station {code}: {e.__class__.__name__}"
                )
                continue

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


def retrieve_waveforms_including_preprocessed(
    station_codes,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    processed_data,
):
    processed_waveforms = processed_data["directory"]
    if processed_waveforms:
        pr_wf_kind = processed_data["wf_kind"]
        pr_wf_factor = processed_data["wf_factor"]
        processed_station_files = processed_data["station_files"]

    retrieved_waveforms = {}

    def handle_station(code):
        """Handle a single station: read preprocessed or retrieve raw waveforms."""
        network, station = code.split(".")

        # Check for preprocessed waveforms
        if processed_waveforms and code in processed_station_files:
            st_obs = read(processed_station_files[code])  # Read preprocessed waveform
            # Convert between waveform kinds if needed
            kind_mapping = {"acceleration": 0, "velocity": 1, "displacement": 2}
            kind_difference = kind_mapping[kind_vd] - kind_mapping[pr_wf_kind]
            operation = (
                st_obs.integrate if kind_difference > 0 else st_obs.differentiate
            )
            for _ in range(abs(kind_difference)):
                operation()
            for tr in st_obs:
                tr.data *= pr_wf_factor  # Apply scaling factor
            return code, st_obs

        # Retrieve raw waveforms if no preprocessed data is available
        network_station_tmp = {network: [station]}
        retrieved_waveforms_tmp = retrieve_waveforms(
            network_station_tmp,
            client_name,
            kind_vd,
            path_observations,
            starttime,
            endtime,
        )
        if code in retrieved_waveforms_tmp:
            return code, retrieved_waveforms_tmp[code]

        # If nothing is retrieved, return None
        return code, None

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=10) as executor:  # Adjust max_workers as needed
        # Submit tasks for all stations
        futures = {
            executor.submit(handle_station, code): code for code in station_codes
        }

        # Collect results as tasks complete
        for future in as_completed(futures):
            code = futures[future]
            try:
                station_code, waveform = future.result()
                if waveform:
                    retrieved_waveforms[station_code] = waveform
            except Exception as e:
                print(f"Error handling {code}: {e}")

    return retrieved_waveforms

