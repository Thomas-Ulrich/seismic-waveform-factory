#!/usr/bin/env python3
import numpy as np
import argparse
import h5py
import matplotlib.pyplot as plt
import matplotlib

parser = argparse.ArgumentParser(
    description="compute Moment Tensor from fault output file"
)
parser.add_argument("filename", help="Moment Tensor (h5)")
parser.add_argument(
    "--idSTF", nargs="+", help="list of the STF to visualize (1st = 0); ", type=int
)
parser.add_argument(
    "--normalized",
    dest="normalized",
    action="store_true",
    help="show normalized stf (integral=1)",
)
parser.add_argument(
    "--cumulative",
    dest="cumulative",
    action="store_true",
    help="show a cumulative plot of all STF",
)
parser.add_argument(
    "--MRF", nargs=1, metavar=("MomentRateFile"), help="save MomentRate to file"
)
parser.add_argument(
    "--time_range",
    nargs=2,
    metavar=("minT", "maxT"),
    help="time range for time slices",
    type=float,
)
parser.add_argument(
    "--extension", nargs=1, default=(["svg"]), help="extension output file"
)

args = parser.parse_args()

fig = plt.figure(figsize=(4.5, 4.5 * 8.0 / 16), dpi=80)
ax = fig.add_subplot(111)

h5f = h5py.File(args.filename, "r")
STF = h5f["normalized_moment_rates"][:, :]
print((STF.shape))
dt = h5f["dt"][0]
nsources, ndt = STF.shape

if not args.normalized:
    moment_tensors = h5f["moment_tensors"][:, :]
    Mom = np.linalg.norm(moment_tensors, axis=1)
    c = np.diag(Mom)
    STF = np.dot(c, STF)
h5f.close()

if not args.idSTF or args.cumulative:
    print("plotting all STF")
    args.idSTF = list(range(0, nsources))

time = np.linspace(0, (ndt - 1) * dt, ndt)
cols = "bgrcykb"


cmap = matplotlib.colormaps["twilight"]
if args.time_range:
    time = time - args.time_range[0]

for ids, i in enumerate(args.idSTF):
    if (args.cumulative) and (i > 0):
        STF[i, :] = STF[i - 1, :] + STF[i, :]
    # plt.plot(time, STF[i,:], cols[ids%7])
    ax.plot(time, STF[i, :], color=cmap(ids / nsources))
if args.time_range:
    ax.set_xlim([0, args.time_range[1] - args.time_range[0]])
else:
    ax.set_xlim([0, ndt * dt])
ax.set_ylim(bottom=0.0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_ylabel("moment rate (Nm/s)")
ax.set_xlabel("time (s)")


prefix = args.filename.split(".")[-2]
fname = f"STF_{prefix}.{args.extension[0]}"
plt.savefig(fname, bbox_inches="tight")
print(f"done writing {fname}")
plt.show()


if args.MRF:
    print(f"saving moment rate function to {args.MRF[0]}")
    if not args.cumulative:
        for i in range(0, nsources):
            STF[i, :] = STF[i - 1, :] + STF[i, :]
    np.savetxt(args.MRF[0], np.column_stack((time, STF[nsources - 1, :])))
