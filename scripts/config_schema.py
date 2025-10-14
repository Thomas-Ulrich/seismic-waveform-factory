CONFIG_SCHEMA = {
    "general": {
        "type": dict,
        "schema": {
            "setup_name": {
                "default": "synthetics",
                "type": str,
                "doc": "Name for the setup, used e.g. in output files.",
            },
            "figure_extension": {
                "default": "png",
                "type": str,
                "doc": "Extension for saved figures.",
            },
            "font_size": {
                "default": 12,
                "type": int,
                "doc": "Font size for plot text.",
            },
            "relative_offset": {
                "default": 0.0,
                "type": float,
                "doc": "Relative vertical offset between traces in plots.",
            },
            "global_legend_labels": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": "Optional list of legend labels for all plots.",
            },
            "projection": {
                "default": None,
                "type": str,
                "doc": "Projection string for coordinate conversion. Optional.",
            },
            "hypocenter": {
                "type": dict,
                "schema": {
                    "lon": {
                        "default": None,
                        "type": float,
                        "doc": "Longitude of the hypocenter.",
                    },
                    "lat": {
                        "default": None,
                        "type": float,
                        "doc": "Latitude of the hypocenter.",
                    },
                    "depth_in_km": {
                        "default": None,
                        "type": float,
                        "doc": "Hypocenter depth in kilometers.",
                    },
                    "onset": {
                        "default": None,
                        "type": str,
                        "doc": "Onset time (formatted as UTC string).",
                    },
                },
            },
            "station_file": {
                "default": None,
                "type": str,
                "doc": "Optional csv file containing station information.",
            },
            "client": {
                "default": "iris",
                "type": str,
                "doc": "obspy client name to retrieve station data",
            },
            "line_colors": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": "List of line colors for plots.",
            },
            "line_widths": {
                "default": [],
                "type": list,
                "element_type": float,
                "doc": "List of line widths for plots.",
            },
            "path_observations": {
                "default": "./observations",
                "type": str,
                "doc": "Path to observed waveform files.",
            },
            "fault_file": {
                "default": "",
                "type": str,
                "doc": (
                    "fault paraview file. used for ploting fault edges, or ",
                    "computing station distance to fault",
                ),
            },
            "misfit": {
                "default": "min_shifted_normalized_rms",
                "type": str,
                "choices": [
                    "min_shifted_normalized_rms",
                    "normalized_rms",
                    "cross-correlation",
                    "time-frequency",
                ],
                "doc": "Type of misfit metric for comparison.",
            },
        },
    },
    "synthetics": {
        "default": [],
        "type": list,
        "doc": (
            "List of synthetic waveform generator configurations. "
            "Each entry must have a 'type' (axitra, instaseis, pyprop8, seissol) "
            "and corresponding parameters."
        ),
        "schema": {
            "type": {
                "default": None,
                "type": str,
                "choices": ["axitra", "instaseis", "pyprop8", "seissol"],
                "doc": "Type of synthetic generator.",
            },
            "name": {
                "default": "",
                "type": str,
                "doc": (
                    "name identifying the synthetics, to be used in ",
                    "waveform_plots.synthetics field",
                ),
            },
            "path": {
                "default": ["."],
                "type": list,
                "element_type": str,
                "doc": "paths where to look for Axitra",
                "when": {"type": "axitra"},
            },
            "fmax": {
                "default": None,
                "type": float,
                "doc": "Maximum frequency.",
                "when": {"type": ("axitra", "pyprop8")},
            },
            "duration": {
                "default": None,
                "type": float,
                "doc": "Duration of synthetics.",
                "when": {"type": ("axitra", "pyprop8")},
            },
            "velocity_model": {
                "default": None,
                "type": str,
                "doc": "Velocity model file.",
                "when": {"type": ("axitra", "pyprop8")},
            },
            "source_files": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": "List of source files or directories containing sources.",
                "when": {"type": ("axitra", "instaseis", "pyprop8")},
            },
            "db": {
                "default": None,
                "type": str,
                "doc": "Database name or path.",
                "when": {"type": "instaseis"},
            },
            "mode": {
                "default": "classical",
                "type": str,
                "choices": ["classical", "green_functions"],
                "doc": "Mode for generating Instaseis synthetics.",
                "when": {"type": "instaseis"},
            },
            "path_computed_synthetics": {
                "default": "./computed_synthetics",
                "type": str,
                "doc": (
                    "Path where to store computed instaseis synthetics ",
                    "and green functions",
                ),
                "when": {"type": "instaseis"},
            },
            "outputs": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": "List of SeisSol prefix_path.",
                "when": {"type": "seissol"},
            },
        },
    },
    "processed_waveforms": {
        "type": dict,
        "schema": {
            "directory": {
                "default": "",
                "type": str,
                "doc": (
                    "Directory containing processed (observed) waveform ",
                    "files, e.g. from Hinet.",
                ),
            },
            "wf_kind": {
                "default": "velocity",
                "type": str,
                "choices": ["acceleration", "velocity", "displacement"],
                "doc": "Kind of processed waveform.",
            },
            "wf_factor": {
                "default": 1.0,
                "type": float,
                "doc": "Scaling factor for processed waveforms.",
            },
        },
    },
    "waveform_plots": {
        "type": list,
        "schema": {
            "type": {
                "default": "p",
                "type": str,
                "choices": ["p", "sh", "generic"],
                "doc": "Type of waveform plot.",
            },
            "enabled": {
                "default": True,
                "type": bool,
                "doc": "Whether this plot is enabled.",
            },
            "misfit": {
                "default": "auto",
                "type": str,
                "choices": [
                    "auto",
                    "min_shifted_normalized_rms",
                    "normalized_rms",
                    "cross-correlation",
                    "time-frequency",
                ],
                "doc": (
                    "Type of misfit metric used. If auto, the misfit defined in ",
                    "general with be used.",
                ),
            },
            "synthetics": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": "list of names of synthetics to be plotted, see synthetics.name",
            },
            "kind": {
                "default": None,
                "type": str,
                "choices": ["acceleration", "velocity", "displacement"],
                "doc": "Kind/type of synthetics to generate (e.g. 'velocity').",
            },
            "stations": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": "List of station codes.",
            },
            "t_before": {
                "default": 0.0,
                "type": (float, dict),
                "doc": (
                    "Time (seconds) before reference time, defining the beginning of ",
                    " the plotted window. Reference time may be user-defined, event ",
                    "origin, or P/S onset depending on plot type. Can be a float or a",
                    " dict of station-specific values (with optional 'default').",
                ),
            },
            "t_after": {
                "default": None,
                "type": (float, dict),
                "doc": (
                    "Time (seconds) after reference time, defining the end of ",
                    " the plotted window. Reference time may be user-defined, event ",
                    "origin, or P/S onset depending on plot type. Can be a float or a",
                    " dict of station-specific values (with optional 'default').",
                ),
            },
            "taper": {"default": True, "type": bool, "doc": "Apply taper to waveform."},
            "normalize": {
                "default": False,
                "type": bool,
                "doc": "Whether to normalize waveform amplitudes.",
            },
            "shift_match_correlation": {
                "default": False,
                "type": bool,
                "doc": "Apply shift to maximize correlation for goodness-of-fit.",
            },
            "scaling": {
                "default": 1.0,
                "type": float,
                "doc": "Scaling factor for waveform amplitudes.",
            },
            "filter_tmin": {
                "default": None,
                "type": float,
                "doc": "Minimum period (seconds) for filter.",
            },
            "filter_tmax": {
                "default": None,
                "type": float,
                "doc": "Maximum period (seconds) for filter.",
            },
            "fault_strike": {
                "default": None,
                "type": float,
                "doc": "Fault strike angle, required for certain FN and FP components.",
            },
            "ncol_per_component": {
                "default": 1,
                "type": int,
                "doc": "Number of plot columns per waveform component.",
            },
            "annotations": {
                "type": dict,
                "schema": {
                    "fields": {
                        "default": ["distance", "azimuth", "misfit"],
                        "type": list,
                        "element_type": str,
                        "doc": "List of annotation fields to show on plots.",
                    },
                    "distance_unit": {
                        "default": "auto",
                        "type": str,
                        "choices": ["auto", "km", "degree"],
                        "doc": "Unit for distance annotations.",
                    },
                },
            },
            "components": {
                "default": [],
                "type": list,
                "element_type": str,
                "doc": (
                    "List of waveform components to plot (e.g. ['Z'], ",
                    "['E', 'N', 'Z'], ['T'], or ['normal', 'parallel'] ",
                    "for fault normal and parallel).",
                ),
            },
        },
    },
}
