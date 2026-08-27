import networkx as nx
import numpy as np
import pyvista as pv
import seissolxdmf
import vtk
from shapely.geometry import Polygon
from vtk.util import numpy_support


def create_vtk_grid(
    xyz: np.ndarray, connect: np.ndarray
) -> vtk.vtkPolyData or vtk.vtkUnstructuredGrid:
    """Create a VTK grid from HDF5 data.

    Parameters:
    xyz (np.ndarray): Node coordinates
    connect (np.ndarray): Connectivity array

    Returns:
    vtk.vtkUnstructuredGrid: The created VTK grid
    """

    n_elements, ndim2 = connect.shape
    grid_type = {3: vtk.vtkPolyData, 4: vtk.vtkUnstructuredGrid}[ndim2]
    grid = grid_type()

    points = vtk.vtkPoints()
    points.SetData(numpy_support.numpy_to_vtk(xyz))
    grid.SetPoints(points)

    cells = vtk.vtkCellArray()
    connect2 = np.zeros((n_elements, ndim2 + 1), dtype=np.int64)
    # number of points in the cell
    connect2[:, 0] = ndim2
    connect2[:, 1:] = connect
    cells.SetCells(n_elements, numpy_support.numpy_to_vtkIdTypeArray(connect2))
    if ndim2 == 3:
        grid.SetPolys(cells)
    else:
        grid.SetCells(vtk.VTK_TETRA, cells)

    return grid


def compute_shapely_polygon(faultfname):
    sx = seissolxdmf.seissolxdmf(faultfname)
    xyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    tags = sx.Read1dData("fault-tag", sx.nElements, isInt=True).astype(np.int32)

    unique_tags = np.unique(tags)
    polygons = []

    for tag in unique_tags:
        # 1. Filter elements matching the current tag
        tag_mask = tags == tag
        tag_connect = connect[tag_mask]

        if tag_connect.size == 0:
            continue

        # 2. Extract only the points used by this specific tag subset
        unique_pt_indices, remapped_connect = np.unique(
            tag_connect, return_inverse=True
        )
        tag_xyz = xyz[unique_pt_indices]
        remapped_connect = remapped_connect.reshape(tag_connect.shape)

        # 3. Create PyVista mesh for the tag component
        grid = create_vtk_grid(tag_xyz, remapped_connect)
        mesh = pv.wrap(grid)

        # 4. Extract boundary edges
        edges = mesh.extract_feature_edges(
            boundary_edges=True,
            feature_edges=False,
            manifold_edges=False,
            non_manifold_edges=False,
        )

        if edges.n_cells == 0:
            continue

        lines = edges.lines.reshape((-1, 3))[:, 1:]

        # 5. Build graph and extract the boundary loop
        G = nx.Graph()
        G.add_edges_from(lines)

        try:
            cycle = nx.find_cycle(G)
        except nx.NetworkXNoCycle:
            continue

        node_indices = [edge[0] for edge in cycle] + [cycle[0][0]]
        poly_points = [edges.points[i, 0:2] for i in node_indices]

        poly = Polygon(poly_points)
        polygons.append(poly)

    return polygons


def get_fault_slip_coords(faultfname):
    sx = seissolxdmf.seissolxdmf(faultfname)
    xyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    xyzcenter = (
        xyz[connect[:, 0], :] + xyz[connect[:, 1], :] + xyz[connect[:, 2], :]
    ) / 3.0
    fault_slip = sx.ReadData("ASl", sx.ReadNdt() - 1)
    # Find elements where fault slip exceeds the threshold
    fault_slip_mask = fault_slip > 0.01
    if not np.any(fault_slip_mask):
        print("No fault slip detected, returning the full coordinate array")
        return xyzcenter
    else:
        return xyzcenter[fault_slip_mask, :]
