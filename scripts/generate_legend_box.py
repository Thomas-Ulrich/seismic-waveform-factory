#!/usr/bin/env python3
import argparse
from itertools import cycle, islice

import matplotlib.pyplot as plt
from config_loader import ConfigLoader
from config_schema import CONFIG_SCHEMA

parser = argparse.ArgumentParser(description=("generate legend box only"))
parser.add_argument("config_file", help="config file describing event and stations")
args = parser.parse_args()

cfg = ConfigLoader(args.config_file, CONFIG_SCHEMA)

colors = cfg["general"]["line_colors"]
line_widths = cfg["general"]["line_widths"]
legends = cfg["general"]["global_legend_labels"]
font_size = cfg["general"]["font_size"]

colors = list(islice(cycle(colors), len(legends)))
line_widths = list(islice(cycle(line_widths), len(legends)))

# Dummy plot data
x = [0, 1]
y = [0, 1]

# Create a figure and a plot
fig, ax = plt.subplots()

# Plot dummy data to generate legend entries
for i, leg in enumerate(legends):
    ax.plot(x, y, color=colors[i], linewidth=line_widths[i], label=leg)

# Add the legend
legend = ax.legend(fontsize=font_size, frameon=False)

# Hide plot elements
ax.set_axis_off()
for line in ax.get_lines():
    line.set_visible(False)

# Save or show the legend box
fname = "legend_box.svg"
plt.savefig(fname, bbox_inches="tight", pad_inches=0.1)
print(f"done generating {fname}")
