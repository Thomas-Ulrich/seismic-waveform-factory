from obspy import read
from lxml.etree import XMLSyntaxError
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.clients.fdsn import Client, RoutingClient
import gzip
import os


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
        max_retries = 5
        retry_count = 0
        while retry_count < max_retries:
            try:
                if retry_count == 0:
                    print(f"getting {level}s for network {network}...")
                inventory = client.get_stations(
                    network=network,
                    station=",".join(stations),
                    level=level,
                    starttime=t1,
                )
                break
            except FDSNNoDataException:
                print(
                    f"FDSNNoDataException Error occurred in get_station for network {network}, with stations {stations}"
                )
                retry_count += 1
                if retry_count == max_retries:
                    print(f"Max retry count reached for {network}. Skipping.")
                    break
            except XMLSyntaxError:
                print(f"XML Syntax Error occurred in get_station for network {network}")
                retry_count += 1
                if retry_count == max_retries:
                    print(f"Max retry count reached for {network}. Skipping.")
                    break
            except gzip.BadGzipFile:
                print(
                    f"gzip.BadGzipFile Error occurred in get_station for network {network}"
                )
                retry_count += 1
                if retry_count == max_retries:
                    print(f"Max retry count reached for {network}. Skipping.")
                    break
        if retry_count == max_retries:
            continue
        for station in inventory[0]:
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

            max_retries = 5
            retry_count = 0
            got_data = False
            while retry_count < max_retries:
                try:
                    print(network, station.code, f"{selected_band}*", t1, t2)
                    if is_routing_client:
                        kwargs = {}
                    else:
                        kwargs = {"attach_response": True}
                    st_obs0 = client.get_waveforms(
                        network=network,
                        station=station.code,
                        location="*",
                        channel=f"{selected_band}*",
                        starttime=t1,
                        endtime=t2,
                        **kwargs,
                    )
                    got_data = True
                    if len(st_obs0) == 0:
                        print(f"got empty stream for {network} {station.code}")
                        got_data = False
                    break
                except FDSNNoDataException:
                    print(f"no waveform available for {station.code}")
                    break
                except XMLSyntaxError:
                    print(f"XML Syntax Error in get_waveforms for {station.code}")
                    retry_count += 1
                    if retry_count == max_retries:
                        print(f"Max retry count for {station.code}. Skipping.")
                        break

            if not got_data:
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
