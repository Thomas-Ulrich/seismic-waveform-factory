

Quick Start
-----------

1. **Installing the package**

    Clone the repository and install the package with:

   .. code-block:: bash

      git clone https://github.com/Thomas-Ulrich/seismic-waveform-factory.git
      cd seismic-waveform-factory
      pip install -e .

   You can now verify that you have access to the CLI by running:

   .. code-block:: bash

      swf --help

   This should output something like:

   .. code-block:: bash

    usage: swf [-h] {compute-multi-cmt,get-focal-mech,plot-map-cmt,plot-stf-cmt,select-stations,plot-waveforms,gen-seissol-sta,gen-map-waveforms,plot-station-map,gen-legend-box} ...

    Seismic Waveform Factory CLI

    options:
      -h, --help            show this help message and exit

    subcommands:
      {compute-multi-cmt,get-focal-mech,plot-map-cmt,plot-stf-cmt,select-stations,plot-waveforms,gen-seissol-sta,gen-map-waveforms,plot-station-map,gen-legend-box}
        compute-multi-cmt   Compute multi-CMT representation from SeisSol fault output file
        get-focal-mech      Extract focal mechanism from multi-CMT file
        plot-map-cmt        Plot a map including beachballs from multi-CMT file
        plot-stf-cmt        Plot source time functions from a multi-CMT file
        select-stations     Select stations ensuring optimal coverage from configuration file
        plot-waveforms      Generate waveform figures from configuration file
        gen-seissol-sta     Generate SeisSol station file from configuration file
        gen-map-waveforms   Generate map with waveforms from configuration file
        plot-station-map    Plot station map from configuration file
        gen-legend-box      Generate legend box from configuration file



2. **CLI usage**


   Example configuration files (`*.yaml`) and usage examples (`run_test.sh`) are provided in the repository, folder `test`.
