import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator
import h5py
import sys


def compute_rigidity(muDescription, muValue, coords):
    """Compute rigidity given model."""
    if muDescription == 0:
        print(f"using constant G: {muValue}")
        return np.ones((coords.shape[0])) * float(muValue)

    elif muDescription == 1:
        print(f"loading 1D velocity file {muValue}")
        depthmu = np.loadtxt(muValue)
        depthmu = depthmu[depthmu[:, 0].argsort()]
        return np.interp(coords[:, 2], depthmu[:, 0], depthmu[:, 1])

    elif muDescription == 2:
        print(f"loading G from 3D netcdf file: {muValue}")
        return getG(muValue, coords)

    elif muDescription == 3:
        print(f"loading 1D velocity file (Sumatra model): {muValue}")
        depthmu = np.loadtxt(muValue)
        depthmu = depthmu[depthmu[:, 0].argsort()]
        y1, y2 = 250000, 1550000
        alpha = np.maximum(0.0, np.minimum(1.0, (coords[:, 1] - y1) / (y2 - y1)))
        depth_seafloor = -5000.0 + 2000.0 * alpha
        print(depth_seafloor)
        return np.interp(coords[:, 2] - depth_seafloor, depthmu[:, 0], depthmu[:, 1])


def compute_slices_array_enforcing_dx(x1, fault_slip, dx, slip_threshold):
    """Compute a linearly spaces coordinate array (for latter slicing the fault output
    along vH and z)"""
    # Useful if large portion of the fault not ruptured
    if np.amax(fault_slip) / slip_threshold < 20:
        print("Warning, slip_threshold is a significant portion of the maximum slip")
    if np.amax(fault_slip) <= slip_threshold:
        print(f"fault slip is below slip_threshold ({slip_threshold:.2f} m) exiting...")
        sys.exit(0)
    x1_slip = x1[fault_slip > slip_threshold]
    x1min = np.amin(x1_slip)
    x1max = np.amax(x1_slip)
    dx1 = x1max - x1min
    nslices = int(dx1 // (dx * 1e3)) + 1
    myrange = np.linspace(0.0, 1.0001, (nslices + 1))
    return np.array([x1min + v0 * dx1 for v0 in myrange])


def compute_slices_array(x1, fault_slip, nslices, slip_threshold):
    """Compute a linearly spaces coordinate array (for latter slicing the fault output
    along vH and z)"""
    # Useful if large portion of the fault not ruptured
    if np.amax(fault_slip) / slip_threshold < 20:
        print("Warning, slip_threshold is a significant portion of the maximum slip")
    x1_slip = x1[fault_slip > slip_threshold]
    x1min = np.amin(x1_slip)
    x1max = np.amax(x1_slip)
    dx1 = x1max - x1min
    myrange = np.linspace(0.0, 1.0001, (nslices + 1))
    return np.array([x1min + v0 * dx1 for v0 in myrange])


def NED2RTP(aMomentTensor):
    """convert array of moment tensor in NED to RTP
    derived from source code of obspy
     https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
    """
    signs = [1, 1, 1, 1, -1, -1]
    indices = [2, 0, 1, 4, 5, 3]
    return np.array(
        [sign * aMomentTensor[:, ind] for sign, ind in zip(signs, indices)]
    ).T


def RTP2NED(aMomentTensor):
    """convert RTP moment tensor to NED
    derived from source code of obspy
    https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
    """
    signs = [1, 1, 1, -1, 1, -1]
    indices = [1, 2, 0, 5, 3, 4]
    return np.array([sign * aMomentTensor[ind] for sign, ind in zip(signs, indices)]).T


def computeMomentTensor(FaceMomentTensor):
    return np.sum(FaceMomentTensor, axis=1)


def compute_seismic_moment(MomentTensor):
    """Compute M0 given the moment tensor (6 components)"""
    fullMomentTensor = np.zeros((3, 3))
    fullMomentTensor[0, :] = [MomentTensor[0], MomentTensor[3], MomentTensor[4]]
    fullMomentTensor[1, :] = [MomentTensor[3], MomentTensor[1], MomentTensor[5]]
    fullMomentTensor[2, :] = [MomentTensor[4], MomentTensor[5], MomentTensor[2]]
    # https://gfzpublic.gfz-potsdam.de/rest/items/item_272892/
    # component/file_541895/content (p6)
    # Note that Mom from the moment rate is unprecise
    # Therefore we compute the moment from the final slip here using the Frobenius norm
    return np.linalg.norm(fullMomentTensor) * np.sqrt(0.5)


def write_point_source_file(fname, point_sources, dt, proj):
    """Write h5 file describing a multi point source model."""
    tags = point_sources.keys()

    for i, fault_tag in enumerate(tags):
        if i == 0:
            MomentTensor = point_sources[fault_tag]["moment_tensors"]
            NormMomentRate = point_sources[fault_tag]["moment_rate_functions"]
            xyz = point_sources[fault_tag]["locations"]
            nsources, ndt = NormMomentRate.shape
            FaultTags = np.zeros((nsources,), dtype=int) + fault_tag
        else:
            MomentTensor = np.concatenate(
                (MomentTensor, point_sources[fault_tag]["moment_tensors"]), axis=0
            )
            NormMomentRate = np.concatenate(
                (NormMomentRate, point_sources[fault_tag]["moment_rate_functions"]),
                axis=0,
            )
            xyz = np.concatenate((xyz, point_sources[fault_tag]["locations"]), axis=0)
            nsources, ndt = point_sources[fault_tag]["moment_rate_functions"].shape
            FaultTags_segment = np.zeros((nsources,), dtype=int) + fault_tag
            FaultTags = np.concatenate((FaultTags, FaultTags_segment), axis=0)

    nsources, ndt = NormMomentRate.shape

    with h5py.File(fname, "w") as h5f:
        nsources, ndt = NormMomentRate.shape
        h5f.create_dataset("NormalizedMomentRate", (nsources, ndt), dtype="d")
        h5f.create_dataset("xyz", (nsources, 3), dtype="d")
        h5f.create_dataset("MomentTensor", (nsources, 6), dtype="d")
        h5f.create_dataset("FaultTags", (nsources,), dtype=int)
        h5f.create_dataset("dt", (1,), dtype="d")
        h5f["dt"][0] = dt
        h5f["NormalizedMomentRate"][:, :] = NormMomentRate[:, :]
        h5f["MomentTensor"][:, :] = MomentTensor[:, :]
        h5f["xyz"][:, :] = xyz[:, :]
        h5f["FaultTags"][:] = FaultTags

        if proj:
            h5f.attrs["CoordinatesConvention"] = np.bytes_("geographic")
        else:
            h5f.attrs["CoordinatesConvention"] = np.bytes_("projected")
    print(f"done writing {fname}")


def getG(asagiFile, xyz):
    """Read 3D Netcdf (Asagi) file and get G at coordinates."""
    fh = Dataset(asagiFile, mode="r")

    if "u" in fh.variables.keys():
        lvar = ["u", "v", "w"]
    elif "x" in fh.variables.keys():
        lvar = ["x", "y", "z"]
    else:
        raise ValueError("wrong variables in netcdf")

    x = fh.variables[lvar[0]][:]
    y = fh.variables[lvar[1]][:]
    z = fh.variables[lvar[2]][:]
    data = fh.variables["data"]
    data = np.swapaxes(data, 0, 2)
    data = [test[1] for test in data.flatten()]
    data = np.array(data)
    data = data.reshape((x.shape[0], y.shape[0], z.shape[0]))
    if z[1] < z[0]:
        data = np.flip(data, axis=2)
        z = np.flip(z)
    fn = RegularGridInterpolator((x, y, z), data)
    return fn(xyz)
