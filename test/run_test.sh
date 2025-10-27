#!/bin/bash

swf compute-multi-cmt spatial data/dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest-fault.xdmf 25e9 --use_geometric_center --slip_threshold " -1e9"  --DH 30 --NZ 1 --proj "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=78.65 +lat_0=41.26"

swf plot-stf-cmt PointSourceFile_dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest_dx30.0_nz1.h5

swf get-focal-mech PointSourceFile_dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest_dx30.0_nz1.h5
swf plot-map-cmt PointSourceFile_dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest_dx30.0_nz1.h5

swf gen-legend-box waveforms_config_regional_sources.yaml

swf select-stations waveforms_config_regional_sources.yaml 14 7
swf select-stations waveforms_config_teleseismic_sources.yaml 8 0 --channel "B*" --station_kind global

swf plot-waveforms waveforms_config_regional_sources.yaml
swf plot-waveforms waveforms_config_teleseismic_sources.yaml

swf gen-seissol-sta waveforms_config_regional_sources.yaml

swf plot-station-map  waveforms_config_regional_sources.yaml
swf plot-station-map waveforms_config_teleseismic_sources.yaml

swf gen-map-waveforms waveforms_config_regional_sources.yaml 0

../src/seismic_waveform_factory/scripts/source_parameters_2_moment_tensor.py

# These 2 latter scripts are not tested yet
#../src/seismic_waveform_factory/scripts/modify_wasp_strong_motion_waves.py
#../src/seismic_waveform_factory/scripts/copy_selected_sac_data_to_folder.py
