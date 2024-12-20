"""
from retrieve_waveforms import retrieve_waveforms
from obspy import UTCDateTime, read
network_station = {'PS': ['INU']}
network_station = {'KIK': ['IGWH']}
client_name = "iris-federator"
kind_vd = "acceleration"
path_observations = "./observations"
starttime= UTCDateTime("2024-01-01T07:08:29.000000Z")
endtime= UTCDateTime("2024-01-01T07:16:49.000000Z")
"""

import argparse
import configparser
from obspy import UTCDateTime
import os
from waveform_figure_utils import compile_station_coords_main
from retrieve_waveforms import retrieve_waveforms_including_preprocessed
import numpy as np

parser = argparse.ArgumentParser(
    description=("prepare obs waveforms to expected LinSlip format")
)
parser.add_argument("config_file", help="config file describing event and stations")
args = parser.parse_args()

config = configparser.ConfigParser()
assert os.path.isfile(args.config_file), f"{args.config_file} not found"
config.read(args.config_file)

setup_name = config.get("GENERAL", "setup_name")
projection = config.get("GENERAL", "projection", fallback=None)

client_name = config.get("GENERAL", "client")
onset = config.get("GENERAL", "onset")
t1 = UTCDateTime(onset)
t_obs_before, t_obs_after = 100, 400
starttime = t1 - t_obs_before
endtime = t1 + t_obs_after

kind_vd = config.get("GENERAL", "kind")

"""
hypo_lon = config.getfloat("GENERAL", "hypo_lon")
hypo_lat = config.getfloat("GENERAL", "hypo_lat")
scaling = config.getfloat("GENERAL", "scaling", fallback=1.0)
normalize = config.getboolean("GENERAL", "normalize", fallback=False)
hypo_depth_in_km = config.getfloat("GENERAL", "hypo_depth_in_km")
"""

station_codes_list = config.get("GENERAL", "stations")
station_codes = [v.strip() for v in station_codes_list.split(",")]

network_station = {}

path_observations = config.get("GENERAL", "path_observations")

processed_data = {}
processed_data["directory"] = config.get(
    "GENERAL", "processed_waveforms", fallback=None
)

if processed_data["directory"]:
    processed_data["wf_kind"] = config.get("GENERAL", "processed_waveforms_kind")
    processed_data["wf_factor"] = config.getfloat(
        "GENERAL", "processed_waveforms_factor", fallback=1.0
    )
    processed_data["station_files"] = get_station_files_dict(
        processed_data["directory"]
    )


station_file = config.get("GENERAL", "station_file", fallback=None)
station_coords = compile_station_coords_main(
    station_codes, station_file, client_name, t1
)

# station_coords = reorder_station_coords_from_azimuth(station_coords, hypo_lon, hypo_lat)

retrieved_waveforms = retrieve_waveforms_including_preprocessed(
    station_codes,
    client_name,
    kind_vd,
    path_observations,
    starttime,
    endtime,
    processed_data,
)
print(retrieved_waveforms)

formatted_time = t1.strftime("%Y %m %d  %H %M %S")
nstation = len(retrieved_waveforms)

out_param = f"""#COMMON ORIGIN (DATE, TIME)
{formatted_time}
#STATIONS
{nstation}
#OUT_DT OUT_SAMPLES ARTIFICIAL_TIMESHIFT
0.4     512         30.
#FILTERS (NUMBER OF FILTERS, PARAMETERS)
2
0.05  0.5
0.00  0.5
--------------------------------------------------------------------------------------------------
LAT       LONG        DATE_BEGIN  TIME_BEGIN    DT     DIVIDE  DETREND FILTER  INTEGRATE  FILENAME
--------------------------------------------------------------------------------------------------
"""

if not os.path.exists("data_ascii_linslip"):
    os.makedirs("data_ascii_linslip")

for code in retrieved_waveforms:
    s_obs = retrieved_waveforms[code]
    trace = s_obs[0]
    dt = trace.stats.delta
    data_to_write = [trace.times(reftime=starttime).T]
    for comp in ["E", "N", "Z"]:
        trace = s_obs.select(component=comp)[0]
        data_to_write += [trace.data]
    data_to_write = np.array(data_to_write).T
    fname = f"data_ascii_linslip/{code}.txt"
    np.savetxt(fname, data_to_write)
    lon, lat = station_coords[code]
    out_param += f"{lon}\t{lat} \t{formatted_time}\t{dt}\t1\t1\t0\t1\t0\t{fname}\n"

print(out_param)
fname = "processseis.in"
with open(fname, "w") as fout:
    fout.write(out_param)
print(f"done writing {fname}")
