#!/usr/bin/env python3
import os
import sys
import importlib.util
from textwrap import indent


def load_schema(schema_path):
    """Dynamically load CONFIG_SCHEMA from a Python file."""
    spec = importlib.util.spec_from_file_location("config_schema", schema_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "CONFIG_SCHEMA"):
        raise ValueError("CONFIG_SCHEMA not found in schema file")
    return module.CONFIG_SCHEMA


def format_type(t):
    if isinstance(t, tuple):
        return ", ".join([format_type(x) for x in t])
    if isinstance(t, type):
        return t.__name__
    return str(t)


def generate_rst_section(name, schema, level=1):
    """Generate RST for one schema section (recursively)."""
    title = f"{name}"
    underline = "-" * len(title) if level == 1 else "~" * len(title)
    rst = [f"{title}\n{underline}\n"]

    if "doc" in schema and schema["doc"]:
        rst.append(schema["doc"] + "\n")

    expected_type = schema.get("type")

    # --- If this is a dict section ---
    if expected_type == dict:
        sub_schema = schema.get("schema", {})
        for key, subschema in sub_schema.items():
            rst.append(generate_rst_entry(key, subschema, level + 1))

    # --- If this is a list section ---
    elif expected_type == list:
        rst.append("\n- Type: list\n")
        if "element_type" in schema:
            rst.append(f"  - Element type: {format_type(schema['element_type'])}\n")
        if "schema" in schema:
            rst.append("- Each element has structure:\n\n")
            for key, subschema in schema["schema"].items():
                rst.append(
                    indent(generate_rst_entry(key, subschema, level + 1), "    ")
                )

    else:
        rst.append(generate_rst_entry(name, schema, level + 1))

    rst.append("\n")
    return "".join(rst)


def generate_rst_entry(key, schema, level=2):
    """Format a single key entry as RST list item."""
    lines = [f"**{key}**\n"]
    t = schema.get("type")
    default = schema.get("default", "None")
    choices = schema.get("choices")
    doc = schema.get("doc", "")

    lines.append(f"  - Type: {format_type(t)}\n")
    lines.append(f"  - Default: {default}\n")
    if choices:
        lines.append(f"  - Choices: {', '.join(map(str, choices))}\n")
    if doc:
        lines.append(f"  - Description: {doc}\n")
    lines.append("\n")
    return "".join(lines)


def main():
    if len(sys.argv) < 3:
        print("Usage: generate_rst_from_schema.py <schema_path> <output_dir>")
        sys.exit(1)

    schema_path = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    schema = load_schema(schema_path)

    doc = [
        "Configuration Parameters\n",
        "=========================\n\n",
        "This document lists all configuration parameters defined in ",
        "CONFIG_SCHEMA.\n\n",
    ]

    for section_name, section_schema in schema.items():
        doc.append(generate_rst_section(section_name, section_schema, level=1))

    output_file = os.path.join(output_dir, "config_parameters.rst")
    with open(output_file, "w") as f:
        f.write("".join(doc))

    print(f"RST documentation written to {output_file}")


if __name__ == "__main__":
    main()
