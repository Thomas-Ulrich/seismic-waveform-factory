# A collection of scripts to generate equivalent multi-cmt solutions, and regional or teleseismic synthetics

1. Generate an equivalent multi-cmt solution with 10 km or less horizontal distance between sources, and  3 rows of layers along deph from fault output file with, e.g.:
```
my_proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=-25.5 +lat_0=-57.5"
prefix=sandwich_10_300_085_resampled
vel_model=SallaresRanero_sf5000.dat
compute_multi_cmt.py spatial --proj "$my_proj" --DH 10 --NZ 3 --slip_threshold 0.4  $prefix-fault.xdmf 1 $vel_model
```
This generate a hdf5 file containing the location, moment tensor and source time function of each sources.
The routines split the faults according to the fault tags, and compute a set of point sources for each segment.

2. The source time function can be plotted with:

```
plot_stf_cmt.py PointSourceFile_${prefix}_dx10.0_nz3.h5
```

Use option ```--beachSize value``` to change the beach ball size.

3. A map of the sources as beach balls can be generated with:

```
plot_map_cmt.py PointSourceFile_${prefix}_10_3.h5 --fault_edge $prefix-fault.xdmf "$my_proj"
```

4. Generate teleseismic synthetics (an example of config file is given in the folder)

```
generate_figure_synthetics.py config.ini
```
