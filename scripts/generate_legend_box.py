#!/usr/bin/env python3
import matplotlib.pyplot as plt
import argparse
import configparser
import os

parser = argparse.ArgumentParser(description=("generate legend box only"))
parser.add_argument("config_file", help="config file describing event and stations")
args = parser.parse_args()


config = configparser.ConfigParser()
assert os.path.isfile(args.config_file), f"{args.config_file} not found"
config.read(args.config_file)

setup_name = config.get("GENERAL", "setup_name")
ext = config.get("GENERAL", "figure_extension")
font_size = config.get("GENERAL", "font_size")

legends = config.get("GENERAL", "legends").split(",")
colors = config.get("GENERAL", "line_colors").split(",")
line_widths = [float(v) for v in config.get("GENERAL", "line_widths").split(",")]
ncolors = len(colors)
if ncolors < len(legends):
    print("enhancing line_colors as not enough colors specified")
    cycol = cycle(colors)
    for i in range(ncolors, len(source_files)):
        colors.append(next(cycol))

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
plt.savefig("legend_box.svg", bbox_inches="tight", pad_inches=0.1)
