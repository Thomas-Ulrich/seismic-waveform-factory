import glob
import gzip
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

from lxml.etree import XMLSyntaxError
from obspy import UTCDateTime, read, read_inventory
from obspy.clients.fdsn import Client, RoutingClient
from obspy.clients.fdsn.header import FDSNException, FDSNNoDataException
from obspy.core.inventory import Inventory
from obspy.core.util.obspy_types import ObsPyException
from requests.exceptions import ConnectionError


def initialize_client(client_name):
    exceptions_to_catch = (UnicodeDecodeError,)
    max_retries = 5

    assert isinstance(client_name, str)

    for retry_count in range(max_retries):
        try:
            if client_name in ["eida-routing", "iris-federator"]:
                return RoutingClient(client_name)
            else:
                return Client(client_name)
        except exceptions_to_catch as e:
            print(f"Error initializing client: {e.__class__.__name__}")
            if retry_count == max_retries - 1:
                raise Exception(
                    f"Max retry count reached while initializing client: {e}"
                )


def load_cached_station_data(network, stations, level, path_observations):
    inv = Inventory()
    stations_not_cached = []
    for station in stations:
        cache_file = os.path.join(path_observations, f"{network}_{station}_{level}.xml")
        if os.path.exists(cache_file):
            print(f"Loading cached data for {network}.{station}")
            station_inv = read_inventory(cache_file)
            inv.extend(station_inv)
        else:
            print(f"No cached data found for {network}.{station}")
            stations_not_cached.append(station)
    return inv, stations_not_cached


def save_station_data(inventory, level, path_observations):
    # Save inventory station-wise regardless of retrieval method
    for net in inventory:
        for sta in net:
            station_inv = Inventory(networks=[net.copy()])
            station_inv[0].stations = [sta]
            station_cache_file = os.path.join(
                path_observations, f"{net.code}_{sta.code}_{level}.xml"
            )
            station_inv.write(station_cache_file, format="STATIONXML")
            print(f"Saved station data {net.code}.{sta.code}")


def filter_channels_by_availability(inventory, starttime, endtime):
    """Filter out channels that don't have data availability in the specified time
    window.

    Parameters:
    -----------
    inventory: obspy.core.inventory.Inventory
        The station inventory to filter
    starttime: str or UTCDateTime
        Start time of the period of interest
    endtime: str or UTCDateTime
        End time of the period of interest

    Returns:
    --------
    obspy.core.inventory.Inventory
        Filtered inventory with only channels that have data in the specified period
    """
    # Convert times to UTCDateTime if they're strings
    if isinstance(starttime, str):
        starttime = UTCDateTime(starttime)
    if isinstance(endtime, str):
        endtime = UTCDateTime(endtime)

    # Create a new inventory
    filtered_networks = []
    code_station_removed = []

    for network in inventory:
        filtered_stations = []

        for station in network:
            filtered_channels = []
            if len(station) == 0:
                print(f"no channel information for {network.code}.{station.code}")
                filtered_stations.append(station)
                continue

            for channel in station:
                # Check if channel has data availability information
                if channel.data_availability:
                    # Check if the requested time window is included
                    # in channel availability
                    chan_start = channel.data_availability.start
                    chan_end = channel.data_availability.end

                    if chan_start <= starttime and chan_end >= endtime:
                        filtered_channels.append(channel)
                else:
                    # keep all channel with no data availability information
                    filtered_channels.append(channel)

            # Only keep station if it has any channels left
            if filtered_channels:
                station.channels = filtered_channels
                filtered_stations.append(station)
            else:
                code_station_removed.append(f"{network.code}.{station.code}")

        # Only keep network if it has any stations left
        if filtered_stations:
            network.stations = filtered_stations
            filtered_networks.append(network)

    # Create new inventory with filtered networks
    filtered_inventory = inventory.copy()
    filtered_inventory.networks = filtered_networks

    if code_station_removed:
        print(
            f"removed {len(code_station_removed)} stations based on data availability"
        )
        print(code_station_removed)
    else:
        station_count_ini = sum(len(network.stations) for network in inventory)
        if station_count_ini > 1:
            print("no station removed based on data availability")

    return filtered_inventory


def get_station_data(
    client_or_clientname,
    network,
    stations,
    channel,
    level,
    t1,
    t2,
    path_observations,
    network_wise=True,
):
    exceptions_to_catch = (
        FDSNException,
        FDSNNoDataException,
        XMLSyntaxError,
        gzip.BadGzipFile,
        ConnectionError,
    )
    os.makedirs(path_observations, exist_ok=True)
    inv, stations_not_cached = load_cached_station_data(
        network, stations, level, path_observations
    )
    if not stations_not_cached:
        return inv

    # Initialize client
    if isinstance(client_or_clientname, str):
        client_or_clientname = initialize_client(client_or_clientname)

    if network_wise:
        station_param = [",".join(stations_not_cached)]
    else:
        station_param = stations_not_cached

    for station in station_param:
        retry_message = f"network {network} at station(s) {station}"
        max_retries = 5
        for retry_count in range(max_retries):
            kargs = {}
            base_url = getattr(client_or_clientname, "base_url", "")
            if "geofon" not in base_url:
                kargs["includeavailability"] = True
            try:
                if retry_count == 0:
                    print(f"Getting {level}s for {retry_message}...")

                inventory = client_or_clientname.get_stations(
                    network=network,
                    station=station,
                    channel=channel,
                    level=level,
                    starttime=t1,
                    endtime=t2,
                    **kargs,
                )

                inv.extend(inventory)
                save_station_data(inventory, level, path_observations)
                break
            except exceptions_to_catch as e:
                if retry_count == max_retries - 1:
                    print(f"Max retry count reached for {retry_message}. Skipping.")
                    print(f"Failed getting {level} for {retry_message}")
                    print(f"Error: {e.__class__.__name__}")
                else:
                    print(
                        f"Error occurred in get_station for {retry_message}:"
                        f" {e.__class__.__name__}"
                    )
                continue
    inv = filter_channels_by_availability(inv, t1, t2)
    return inv


def get_waveforms(client, network, station, selected_band, t1, t2):
    max_retries = 1
    exceptions_to_catch = (
        FDSNException,
        FDSNNoDataException,
        XMLSyntaxError,
        gzip.BadGzipFile,
        ConnectionError,
    )
    st_obs0 = False

    for retry_count in range(max_retries):
        try:
            st_obs0 = client.get_waveforms(
                network=network,
                station=station,
                location="00,",
                channel=f"{selected_band}",
                starttime=t1,
                endtime=t2,
            )
            if not st_obs0:
                print(f"Got empty stream for {network} {station}")
            break
        except exceptions_to_catch as e:
            if retry_count == max_retries - 1:
                print(f"Max retry count reached for {station}. Skipping.")
                print(f"Error: {e.__class__.__name__}")
            else:
                print(
                    f"Error occurred in get_waveforms for {station}: "
                    f"{e.__class__.__name__}"
                )
            continue

    if not st_obs0:
        print(f"No waveform available for {network}.{station}")

    return st_obs0


def get_pre_filt(is_regional, fs):
    nyq = fs / 2
    f1 = 0.001 if is_regional else 0.00033
    f2 = 0.005 if is_regional else 0.001
    f3 = min(15, nyq * 0.8)  # passband end
    f4 = min(30, nyq * 0.95)  # high cut
    return [f1, f2, f3, f4]


def select_band_with_data(
    stream, channels, is_regional, priorities=["B", "H", "E", "M"]
):
    """Select a band based on priority, considering only channels with actual data.

    :param stream: ObsPy Stream object containing retrieved waveform data.
    :param channels: List of channel codes from the inventory.
    :return: The selected band, or None if no suitable band has data.
    """
    instrument_codes = ["N", "H"] if is_regional else ["H"]

    for sub_band in instrument_codes:
        for band in priorities:
            full_band = band + sub_band
            if any(channel.startswith(full_band) for channel in channels):
                # Check if the stream contains data for this band
                if stream.select(channel=f"{full_band}*"):
                    return full_band
    return None


def _retrieve_waveforms(
    network_station,
    client_name,
    kind_vd,
    path_observations,
    t1,
    t2,
    output_format,
    selected_channels,
    is_regional,
):
    level = "response"
    client_initialized = False

    os.makedirs(path_observations, exist_ok=True)
    retrieved_waveforms = {}
    for network, stations in network_station.items():
        sac_file_needed = False
        if output_format == "sac":
            sac_file_needed = True
            for station in stations.copy():
                # Get SAC data files
                pattern = f"{network}.{station}.*.sac"
                search_path = os.path.join(path_observations, pattern)
                matching_files = glob.glob(search_path)

                # Get PZ files
                pattern = f"SAC_PZs_{network}_{station}_*"
                search_path = os.path.join(path_observations, pattern)
                matching_pz_files = glob.glob(search_path)

                if matching_files and matching_pz_files:
                    sac_file_needed = False

        for station in stations.copy():
            pattern = f"{network}.{station}_*_{kind_vd}_{t1.date}.mseed"
            search_path = os.path.join(path_observations, pattern)
            matching_files = sorted(glob.glob(search_path))
            if matching_files:
                fullfname = matching_files[0]
                if not sac_file_needed:
                    print(f"reading the data from {fullfname}")
                    code = f"{network}.{station}"
                    retrieved_waveforms[code] = read(fullfname)
                    stations.remove(station)
                else:
                    print(
                        f"{fullfname} found, but we will redownload the data because "
                        "the sac file is needed (and we need raw data, not instrumented"
                        "  corrected)"
                    )

        if not stations:
            continue
        if not client_initialized:
            client = initialize_client(client_name)

        if level == "channel":
            inventory = get_station_data(
                client,
                network,
                stations,
                selected_channels,
                level,
                t1,
                t2,
                path_observations,
                network_wise=True,
            )
        else:
            inventory = get_station_data(
                client,
                network,
                stations,
                selected_channels,
                level,
                t1,
                t2,
                path_observations,
                network_wise=False,
            )
        if len(inventory) == 0:
            print(f"could not get {level} for {network}")
            continue

        for station in [sta for net in inventory for sta in net]:
            code = f"{network}.{station.code}"
            print(f"requesting waveform data for {network}.{station.code}")
            channels = [channel.code for channel in station]
            print("available channels", channels)
            kargs = {
                "client": client,
                "network": network,
                "station": station.code,
                "selected_band": selected_channels,
                "t1": t1,
                "t2": t2,
            }
            print(kargs)
            # Get waveforms for all channels
            st_obs0 = get_waveforms(**kargs)
            if not st_obs0:
                continue

            st_obs0.attach_response(inventory)

            # Dynamically select a band with actual data
            selected_band = select_band_with_data(st_obs0, channels, is_regional)
            if not selected_band:
                print(f"No valid waveform data at station {station.code}")
                continue

            # Filter the stream for the selected band
            st_obs0 = st_obs0.select(channel=f"{selected_band}*")
            if output_format == "sac":
                # Save as SAC without response removal
                write_sac_files(st_obs0, inventory, path_observations)

            # define a filter band to prevent amplifying noise during the deconvolution
            # pre_filt = [0.00033, 0.001, 1.0, 3.0]
            # Process and save as MSEED with response removal
            output_dic = {
                "acceleration": "ACC",
                "velocity": "VEL",
                "displacement": "DISP",
            }
            pre_filt = get_pre_filt(is_regional, fs=st_obs0[0].stats.sampling_rate)

            def get_water_level(output, instrument_code):
                """
                adapt water_level to instrument measured quantity
                (see warning in)
                https://docs.obspy.org/master/packages/autogen/
                obspy.core.trace.Trace.remove_response.html
                """
                if output == "ACC" and instrument_code == "H":
                    return None
                else:
                    return 60.0

            water_level = get_water_level(
                output_dic[kind_vd], st_obs0[0].stats.channel[1]
            )

            try:
                st_obs0.remove_response(
                    output=output_dic[kind_vd],
                    pre_filt=pre_filt,
                    water_level=water_level,
                    zero_mean=True,
                    taper=True,
                    taper_fraction=0.05,
                    inventory=inventory,
                )
            # except (ValueError, ObsPyException) as e:
            except ObsPyException as e:
                # obspy.core.util.obspy_types.ObsPyException:
                # Can not use evalresp on response with no response stages.
                print(
                    f"Error in st_obs0.remove_response at station {code}:"
                    f" {e.__class__.__name__}"
                )
                continue
            try:
                st_obs0.rotate(method="->ZNE", inventory=inventory)
            except ValueError as e:
                # get rid of this rare error:
                # raise ValueError("The given directions are not linearly independent, "
                # ValueError: The given directions are not linearly independent,
                # at least within numerical precision. Determinant of the base change
                # matrix: 0
                print(
                    f"Error in st_obs0.rotate at station {code}: {e.__class__.__name__}"
                )
                continue

            write_mseed_files(
                st_obs0, code, selected_band, kind_vd, t1, path_observations
            )

            retrieved_waveforms[code] = st_obs0
    return retrieved_waveforms


def write_mseed_files(st_obs0, code, selected_band, kind_vd, t1, path_observations):
    """Write seismic data in MSEED format.

    Parameters:
    -----------
    st_obs0 : obspy.Stream
        Stream object containing the seismic data
    code : str
        Station code
    selected_band : str
        Selected frequency band
    kind_vd : str
        Kind of data (velocity, displacement, etc.)
    t1 : datetime
        Start time
    path_observations : str
        Path where the MSEED files will be saved
    """
    fname = f"{code}_{selected_band}_{kind_vd}_{t1.date}.mseed"
    fullfname = os.path.join(path_observations, fname)
    st_obs0.write(fullfname, format="MSEED")
    print(f"done writing {fullfname}")


def write_sac_files(st_obs0, inventory, path_observations):
    """Write seismic data in SAC format and corresponding SAC pole-zero files.

    Parameters:
    -----------
    st_obs0 : obspy.Stream
        Stream object containing the seismic data
    inventory : obspy.Inventory
        Station inventory containing metadata
    path_observations : str
        Path where the SAC files will be saved
    """
    for tr in st_obs0:
        # Get station coordinates from inventory
        inv_channel = inventory.select(
            network=tr.stats.network,
            station=tr.stats.station,
            channel=tr.stats.channel,
            location=tr.stats.location,
        )

        # Get coordinates from the first (and usually only) selected station
        if inv_channel[0].stations:
            station = inv_channel[0].stations[0]
            channel = inv_channel[0][0][0]

            # Initialize SAC dictionary if it doesn't exist
            if not hasattr(tr.stats, "sac"):
                tr.stats.sac = {}

            # Set station coordinates in SAC header
            tr.stats.sac.stla = station.latitude
            tr.stats.sac.stlo = station.longitude
            tr.stats.sac.stel = station.elevation

            # Set component orientation based on channel code
            comp = tr.stats.channel[-1].upper()
            if comp in ["1", "2"]:
                if hasattr(channel, "azimuth") and hasattr(channel, "dip"):
                    tr.stats.sac.cmpaz = channel.azimuth  # Azimuth from north (degrees)
                    tr.stats.sac.cmpinc = (
                        90 + channel.dip
                    )  # Convert dip to incidence angle
                else:
                    raise ValueError(
                        f"channel.azimuth abd channel.dip unknown for {tr.id}"
                    )

        # Write SAC waveform file
        fname = f"{tr.stats.network}.{tr.stats.station}.{tr.stats.channel}.sac"
        fullfname = os.path.join(path_observations, fname)
        tr.write(fullfname, format="SAC")

        # Write SAC pole-zero file
        location_padded = tr.stats.location.ljust(2, "_")

        pz_fname = (
            f"SAC_PZs_{tr.stats.network}_{tr.stats.station}_"
            f"{tr.stats.channel}_{location_padded}"
        )
        fullfname_pz = os.path.join(path_observations, pz_fname)
        inv_channel.write(fullfname_pz, format="SACPZ")
        print(f"done writing {fullfname_pz}")


def retrieve_waveforms_preprocessed(
    code,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    processed_data,
):
    keys_to_check = {"directory", "wf_kind", "wf_factor", "station_files"}
    if keys_to_check.issubset(processed_data.keys()):
        pr_wf_kind = processed_data["wf_kind"]
        pr_wf_factor = processed_data["wf_factor"]
        processed_station_files = processed_data["station_files"]

        # Check for preprocessed waveforms
        if code in processed_station_files:
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
            return st_obs
    elif processed_data:
        raise ValueError(
            f"{processed_data} is given but does not provides all keys_to_check"
        )
    else:
        return None


def handle_station(
    code,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    output_format,
    processed_data,
    channel,
    is_regional,
):
    """Handle a single station: read preprocessed or retrieve raw waveforms."""
    network, station = code.split(".")

    st_obs = retrieve_waveforms_preprocessed(
        code,
        client_name,
        kind_vd,
        path_observations,
        starttime,
        endtime,
        processed_data,
    )
    if st_obs:
        return code, st_obs

    # Retrieve raw waveforms if no preprocessed data is available
    network_station_tmp = {network: [station]}
    retrieved_waveforms_tmp = _retrieve_waveforms(
        network_station_tmp,
        client_name,
        kind_vd,
        path_observations,
        starttime,
        endtime,
        output_format,
        channel,
        is_regional,
    )
    if code in retrieved_waveforms_tmp:
        return code, retrieved_waveforms_tmp[code]

    # If nothing is retrieved, return None
    return code, None


def retrieve_waveforms(
    station_codes,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    processed_data={},
    output_format="mseed",
    parallel=False,
    channel="*",
    is_regional=True,
):
    retrieved_waveforms = {}
    assert output_format in ["sac", "mseed"]
    # Predefine arguments using partial
    handle_station_partial = partial(
        handle_station,
        client_name=client_name,
        kind_vd=kind_vd,
        path_observations=path_observations,
        starttime=starttime,
        endtime=endtime,
        processed_data=processed_data,
        output_format=output_format,
        channel=channel,
        is_regional=is_regional,
    )

    if parallel:
        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=10) as executor:
            # Submit tasks for all stations
            futures = {
                executor.submit(handle_station_partial, code): code
                for code in station_codes
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

    else:
        for code in station_codes:
            station_code, waveform = handle_station_partial(code)
            if waveform:
                retrieved_waveforms[station_code] = waveform

    return retrieved_waveforms
