from copy import deepcopy

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedSeq


def determine_config_scale(cfg):
    """
    Categorize input yaml as 'global' or 'regional'
    depending on synthetic type.
    """
    is_teleseismic = {
        syn["name"]: (syn["type"] == "instaseis") for syn in cfg["synthetics"]
    }

    teleseismic = False
    regional = False
    for wf_plot in cfg["waveform_plots"]:
        if not wf_plot.get("enabled", True):
            continue
        global_sta = any(is_teleseismic[name] for name in wf_plot["synthetics"])
        if global_sta:
            teleseismic = True
        else:
            regional = True
    return {"regional": regional, "global": teleseismic}


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
        if syn["type"] == "seissol":
            # todo remove hardocded
            durations += [200.0]
        elif syn["type"] != "instaseis":
            durations += [syn["duration"]]
    return durations


def yaml_dump(cfg_ini, fname):
    cfg = deepcopy(cfg_ini)

    yaml = YAML()
    yaml.default_flow_style = None  # use flow only when we say so
    yaml.indent(mapping=2, sequence=4, offset=2)

    def remove_none(d):
        """Recursively remove keys with None values and items that are None."""
        if isinstance(d, dict):
            return {k: remove_none(v) for k, v in d.items() if v is not None}
        elif isinstance(d, list):
            return [remove_none(v) for v in d if v is not None]
        else:
            return d

    def set_flow_for_leaf_lists(node):
        """Recursively mark leaf lists (only scalars) for flow style."""
        if isinstance(node, dict):
            for v in node.values():
                set_flow_for_leaf_lists(v)
        elif isinstance(node, list):
            if all(not isinstance(x, (list, dict)) for x in node):
                # leaf list → make inline
                seq = CommentedSeq(node)
                seq.fa.set_flow_style()
                return seq
            else:
                for i, v in enumerate(node):
                    node[i] = set_flow_for_leaf_lists(v)
        return node

    cfg = remove_none(cfg)
    cfg = set_flow_for_leaf_lists(cfg)

    with open(fname, "w") as f:
        yaml.dump(cfg, f)
    print(f"done creating {fname}")
