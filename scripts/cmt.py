import numpy as np
import h5py
import sys


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


def NED2RTP(moment_tensor):
    """convert array of moment tensor in NED to RTP
    derived from source code of obspy
     https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
    """
    signs = [1, 1, 1, 1, -1, -1]
    indices = [2, 0, 1, 4, 5, 3]
    return np.array(
        [sign * moment_tensor[:, ind] for sign, ind in zip(signs, indices)]
    ).T


def RTP2NED(moment_tensor):
    """convert RTP moment tensor to NED
    derived from source code of obspy
    https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
    """
    signs = [1, 1, 1, -1, 1, -1]
    indices = [1, 2, 0, 5, 3, 4]
    return np.array([sign * moment_tensor[ind] for sign, ind in zip(signs, indices)]).T


def compute_moment_tensor(face_moment_tensor):
    return np.sum(face_moment_tensor, axis=1)


def compute_seismic_moment(moment_tensor):
    """Compute M0 given the moment tensor (6 components)"""
    full_moment_tensor = np.zeros((3, 3))
    full_moment_tensor[0, :] = [moment_tensor[0], moment_tensor[3], moment_tensor[4]]
    full_moment_tensor[1, :] = [moment_tensor[3], moment_tensor[1], moment_tensor[5]]
    full_moment_tensor[2, :] = [moment_tensor[4], moment_tensor[5], moment_tensor[2]]
    # https://gfzpublic.gfz-potsdam.de/rest/items/item_272892/
    # component/file_541895/content (p6)
    # Note that Mom from the moment rate is unprecise
    # Therefore we compute the moment from the final slip here using the Frobenius norm
    return np.linalg.norm(full_moment_tensor) * np.sqrt(0.5)


def write_point_source_file(fname, point_sources, dt, proj, is_potency):
    """Write h5 file describing a multi point source model."""

    if not point_sources:
        raise ValueError("point_sources cannot be empty")

    # Extract all data using list comprehensions
    fault_data_list = list(point_sources.values())

    moment_tensor = np.concatenate(
        [data["moment_tensors"] for data in fault_data_list], axis=0
    )
    normalized_moment_rate = np.concatenate(
        [data["moment_rate_functions"] for data in fault_data_list], axis=0
    )
    xyz = np.concatenate([data["locations"] for data in fault_data_list], axis=0)
    segment_indices = np.concatenate(
        [data["segment_indices"] for data in fault_data_list], axis=0
    )

    fault_tags = np.concatenate(
        [
            np.full(data["moment_rate_functions"].shape[0], tag, dtype=int)
            for tag, data in point_sources.items()
        ],
        axis=0,
    )

    with h5py.File(fname, "w") as h5f:
        h5f.create_dataset("normalized_moment_rate", data=normalized_moment_rate, dtype="d")
        h5f.create_dataset("xyz", data=xyz, dtype="d")
        h5f.create_dataset("moment_tensor", data=moment_tensor, dtype="d")
        h5f.create_dataset("fault_tags", data=fault_tags, dtype=int)
        h5f.create_dataset("segment_indices", data=segment_indices, dtype=int)
        h5f.create_dataset("dt", data=dt, dtype="d")

        convention = "geographic" if proj else "projected"
        h5f.attrs["coordinates_convention"] = convention.encode("utf-8")
        h5f.attrs["is_potency"] = is_potency

    print(f"done writing {fname}")
