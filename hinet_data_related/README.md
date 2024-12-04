# scripts to download and process KIKnet data

```
# download data from hinet, e.g. with 
python get_waveform_event.py --hypocenter 137.27 37.48 --radius_range 1 2 --logging_data user password  --origin_time_utc "2024-01-01 07:10:09
# convert to miniseed
python cnt_to_obspy.py BOSAI_2024-01-01_07:10:09

```
Then add to the ini file:
```
processed_waveforms= path to BOSAI_2024-01-01_07:10:09_miniseed
station_file= path to BOSAI_2024-01-01_07:10:09_miniseed/stations.csv
processed_waveforms_kind = velocity
```
