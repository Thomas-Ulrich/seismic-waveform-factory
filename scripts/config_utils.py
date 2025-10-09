def _categorize_by_scale(cfg, field):
    """
    Helper to categorize waveform data as 'global' or 'regional'
    depending on synthetic type.
    """
    is_teleseismic = {
        syn["name"]: (syn["type"] == "instaseis") for syn in cfg["synthetics"]
    }

    categorized = {"global": set(), "regional": set()}

    for wf_plot in cfg["waveform_plots"]:
        if not wf_plot.get("enabled", True):
            continue
        global_sta = any(is_teleseismic[name] for name in wf_plot["synthetics"])
        key = "global" if global_sta else "regional"

        # wf_plot[field] can be list (stations) or single value (kind)
        values = (
            wf_plot[field] if isinstance(wf_plot[field], list) else [wf_plot[field]]
        )
        categorized[key].update(values)

    # Convert sets to lists for final output
    return {k: list(v) for k, v in categorized.items()}


def categorize_stations_by_scale(cfg):
    """Categorize stations as 'global' or 'regional' based on synthetic model scale."""
    return _categorize_by_scale(cfg, field="stations")


def categorize_waveform_kind_by_scale(cfg):
    """
    Categorize waveform kinds (velocity, acceleration, etc.)
    as 'global' or 'regional'.
    """
    return _categorize_by_scale(cfg, field="kind")


def extract_instaseis_db(cfg):
    """extract instaseis db from config file"""
    dbs = []
    for syn in cfg["synthetics"]:
        if syn["type"] == "instaseis":
            dbs += [syn["db"]]
    return dbs


def extract_regional_durations(cfg):
    """extract instaseis db from config file"""
    durations = []
    for syn in cfg["synthetics"]:
        if syn["type"] != "instaseis":
            durations += [syn["duration"]]
    return durations
