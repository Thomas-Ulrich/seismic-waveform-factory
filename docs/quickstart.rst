

Quick Start
-----------

1. **Generate an equivalent multi-CMT solution**

   .. code-block:: bash

      my_proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=-25.5 +lat_0=-57.5"
      prefix=sandwich_10_300_085_resampled
      vel_model=SallaresRanero_sf5000.dat
      compute_multi_cmt.py spatial --proj "$my_proj" --DH 10 --NZ 3 --slip_threshold 0.4  $prefix-fault.xdmf 1 $vel_model

   This generates an HDF5 file containing location, moment tensor, and source time function for each source. Faults are split by tags and a set of point sources is computed for each segment.

2. **Plot the source time functions of the multi-CMT solution**

   .. code-block:: bash

      plot_stf_cmt.py PointSourceFile_${prefix}_dx10.0_nz3.h5

   Use ``--beachSize value`` to adjust the beach ball size.

3. **Generate a map of the multi-CMT solution sources as beach balls**

   .. code-block:: bash

      plot_map_cmt.py PointSourceFile_${prefix}_10_3.h5 --fault_edge $prefix-fault.xdmf "$my_proj"

4. **Generate waveforms plots**

   .. code-block:: bash

      generate_figure_synthetics.py config.ini

   Example configuration files are provided in the repository.
