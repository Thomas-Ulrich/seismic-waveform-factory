#!/bin/bash

../src/seismic_waveform_factory/compute_multi_cmt.py spatial data/dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest-fault.xdmf 25e9 --use_geometric_center --slip_threshold " -1e9"  --DH 30 --NZ 1 --proj "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=78.65 +lat_0=41.26"

../src/seismic_waveform_factory/plot_stf_cmt.py PointSourceFile_dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest_dx30.0_nz1.h5
../src/seismic_waveform_factory/get_focal_mechanism_multi_cmt.py PointSourceFile_dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest_dx30.0_nz1.h5
../src/seismic_waveform_factory/plot_map_cmt.py PointSourceFile_dyn_0042_coh0.25_1.0_B1.0_C0.3_R0.55_extractedtest_dx30.0_nz1.h5

../src/seismic_waveform_factory/generate_legend_box.py waveforms_config_regional_sources.yaml


../src/seismic_waveform_factory/select_stations.py waveforms_config_regional_sources.yaml 14 7
../src/seismic_waveform_factory/select_stations.py waveforms_config_teleseismic_sources.yaml 8 0 --channel "B*" --station_kind global

../src/seismic_waveform_factory/generate_figure_synthetics.py waveforms_config_regional_sources.yaml
../src/seismic_waveform_factory/generate_figure_synthetics.py waveforms_config_teleseismic_sources.yaml

../src/seismic_waveform_factory/generate_seissol_station_file.py waveforms_config_regional_sources.yaml

../src/seismic_waveform_factory/plot_station_map.py waveforms_config_regional_sources.yaml
../src/seismic_waveform_factory/plot_station_map.py waveforms_config_teleseismic_sources.yaml

../src/seismic_waveform_factory/generate_map_with_waveforms.py waveforms_config_regional_sources.yaml 0

../src/seismic_waveform_factory/source_parameters_2_moment_tensor.py

# These 2 latter scripts are not tested yet
#../src/seismic_waveform_factory/modify_wasp_strong_motion_waves.py
#../src/seismic_waveform_factory/copy_selected_sac_data_to_folder.py
