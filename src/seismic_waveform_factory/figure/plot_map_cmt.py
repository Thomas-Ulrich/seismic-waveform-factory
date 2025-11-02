#!/usr/bin/env python3
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from obspy.imaging.mopad_wrapper import beach
from pyproj import Transformer

from seismic_waveform_factory.utils import cmt
from seismic_waveform_factory.utils.scalebar import scale_bar


def plot_outer_boundary(vertex, connect, ax):
    print("plotting mesh outter boundary")
    import trimesh

    mesh = trimesh.Trimesh(vertex, connect)
    # list vertex of the face boundary
    unique_edges = mesh.edges[
        trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
    ]
    for edge in unique_edges:
        j, k = edge
        xy = mesh.vertices[[j, k], 0:2]
        ax.plot(xy[:, 0], xy[:, 1], "k", linewidth=0.5)


def setup_map(ax, MapBoundaries, grid_spacing, gridlines_left=True):
    """Setup the map with cartopy."""
    ax.set_extent(MapBoundaries, crs=ccrs.PlateCarree())
    scale = "10m"
    ax.add_feature(
        cfeature.LAND.with_scale(scale), facecolor="whitesmoke", rasterized=True
    )
    ax.add_feature(cfeature.OCEAN.with_scale(scale), rasterized=True)
    ax.add_feature(cfeature.COASTLINE.with_scale(scale))
    ax.add_feature(cfeature.BORDERS.with_scale(scale), linestyle=":")
    locs = np.arange(-180, 180, grid_spacing)
    gl = ax.gridlines(draw_labels=True, ylocs=locs, xlocs=locs)
    gl.right_labels = False
    gl.top_labels = False
    gl.left_labels = gridlines_left


def main(args):
    h5f = h5py.File(args.filename, "r")
    moment_tensors = h5f["moment_tensors"][:, :]
    segment_indices = h5f["segment_indices"][:, :]
    fault_tags = h5f["fault_tags"][:]

    xyz = h5f["xyz"][:, :]
    h5f.close()
    xyz[:, 0] += args.x0y0proj[0]
    xyz[:, 1] += args.x0y0proj[1]

    if args.proj:
        # epsg:4226 is geocentric (lat, lon)
        transformer = Transformer.from_crs(args.proj[0], "epsg:4326", always_xy=True)
        print("projecting back to WGS84")
        xyz[:, 0], xyz[:, 1] = transformer.transform(xyz[:, 0], xyz[:, 1])
        print((xyz[:, :]))

    if not args.MapBoundaries:
        print("MapBoundaries not specified, inferring...")
        xmin = np.amin(xyz[:, 0])
        xmax = np.amax(xyz[:, 0])
        ymin = np.amin(xyz[:, 1])
        ymax = np.amax(xyz[:, 1])
        dx = max(ymax - ymin, xmax - xmin, 0.1)
        xmin = xmin - 0.3 * dx
        xmax = xmax + 0.3 * dx
        ymin = ymin - 0.3 * dx
        ymax = ymax + 0.3 * dx
        args.MapBoundaries = [xmin, xmax, ymin, ymax]
    print(
        f"lon_min, lon_max {args.MapBoundaries[0:2]}, \
    lat_min, lat_max {args.MapBoundaries[2:4]}"
    )

    dydx = (args.MapBoundaries[3] - args.MapBoundaries[2]) / (
        args.MapBoundaries[1] - args.MapBoundaries[0]
    )

    fig = plt.figure(figsize=(8, 8 * dydx), dpi=80)
    ax = plt.axes(projection=ccrs.PlateCarree())
    setup_map(ax, args.MapBoundaries, args.dGrid[0])

    if args.fault_edge:
        import seissolxdmf as sx

        sx0 = sx.seissolxdmf(args.fault_edge[0])
        xyz0 = sx0.ReadGeometry()
        connect0 = sx0.ReadConnect()

        transformer = Transformer.from_crs(
            args.fault_edge[1], "epsg:4326", always_xy=True
        )
        print("projecting fault nodes back to WGS84")
        xyz0[:, 0], xyz0[:, 1] = transformer.transform(xyz0[:, 0], xyz0[:, 1])
        plot_outer_boundary(xyz0, connect0, ax)

    cmap = matplotlib.colormaps["twilight"]

    nsrc = moment_tensors.shape[0]
    aMw = np.zeros((nsrc,))
    for isrc in range(nsrc):
        M0all = cmt.compute_seismic_moment(moment_tensors[isrc, :])
        aMw[isrc] = 2.0 / 3.0 * np.log10(M0all) - 6.07

    unique_tags = np.unique(fault_tags)
    for tag in unique_tags:
        print(f"processing segment tagged {tag}...")
        ids = np.where(fault_tags == tag)

        moment_tensors_selected = moment_tensors[ids]
        segment_indices_selected = segment_indices[ids]
        xyz_selected = xyz[ids]
        aMw_seleted = aMw[ids]
        mw_average = np.average(aMw_seleted)

        nsrc = moment_tensors.shape[0]

        for isrc in range(nsrc):
            Mw = aMw_seleted[isrc]

            print(f"{isrc:5d}: {Mw:.2f}", end=" ")
            if isrc % 10 == 9:
                print()
            moment_tensor = cmt.RTP2NED(moment_tensors_selected[isrc, :] / M0all)

            color = cmap(isrc / nsrc)
            if args.unicolor:
                color = "b"

            si = segment_indices_selected[isrc]
            xy = xyz_selected[isrc, 0:2] + np.array([0, args.shift * si[1]])

            beach1 = beach(
                moment_tensor,
                xy=xy,
                width=(Mw / mw_average) ** 2 * args.beachSize[0] * 0.04,
                linewidth=0.2,
                mopad_basis="NED",
                facecolor=color,
            )
            ax.add_collection(beach1)

    if not args.unicolor:
        norm = matplotlib.colors.Normalize(vmin=0, vmax=nsrc)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=0.5 / max(1, dydx), label="id source")

    scale_bar(ax, (0.1, 0.1), args.scalebarSize[0])

    prefix = args.filename.split(".")[-2]
    fname = f"BeachBallPlot_{prefix}.{args.extension[0]}"
    plt.savefig(fname, bbox_inches="tight", dpi=100)
    print(f"done writing {fname}")
