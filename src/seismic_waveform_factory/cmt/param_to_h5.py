import numpy as np
from seismic_waveform_factory.utils import cmt
import argparse
import os

# some code adapted from https://instaseis.net/_modules/instaseis/source.html


def compute_moment_tensor_NED(strike, dip, rake, M0):
    cs = np.cos(strike)
    c2s = np.cos(2.0 * strike)
    cd = np.cos(dip)
    c2d = np.cos(2.0 * dip)
    cl = np.cos(rake)

    ss = np.sin(strike)
    s2s = np.sin(2.0 * strike)
    sd = np.sin(dip)
    s2d = np.sin(2.0 * dip)
    sl = np.sin(rake)

    moment_tensor = np.zeros((6,))
    # 0   1  2  3  4  5
    # xx,yy,zz,xy,xz,yz
    # http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:65580/component/
    # escidoc:65579/IS_3.8_rev1.pdf (eq 5)
    # with x y z : NED
    moment_tensor[0] = -M0 * (sd * cl * s2s + s2d * sl * np.power(ss, 2))
    moment_tensor[1] = M0 * (sd * cl * s2s - s2d * sl * np.power(cs, 2))
    moment_tensor[2] = M0 * (s2d * sl)
    moment_tensor[3] = M0 * (sd * cl * c2s + 0.5 * s2d * sl * s2s)
    moment_tensor[4] = -M0 * (cd * cl * cs + c2d * sl * ss)
    moment_tensor[5] = -M0 * (cd * cl * ss - c2d * sl * cs)
    return moment_tensor


def asymmetric_cosine(trise, tfall=None, npts=10000, dt=0.1):
    """Initialize a source time function with asymmetric cosine, normalized to 1.

    :param trise: rise time
    :type trise: float
    :param tfall: fall time, defaults to trise
    :type trise: float, optional
    :param npts: number of samples
    :type npts: int, optional
    :param dt: sample interval
    :type dt: float, optional
    """
    # initialize
    if not tfall:
        tfall = trise
    t = np.linspace(0, npts * dt, npts, endpoint=False)
    asc = np.zeros(npts)

    # build slices
    slrise = t <= trise
    slfall = np.logical_and(t > trise, t <= trise + tfall)

    # compute stf
    asc[slrise] = 1.0 - np.cos(np.pi * t[slrise] / trise)
    asc[slfall] = 1.0 - np.cos(np.pi * (t[slfall] - trise + tfall) / tfall)

    # normalize
    asc /= trise + tfall

    return asc


parser = argparse.ArgumentParser(
    description="convert file describing kinematic model in param format to h5"
)
parser.add_argument("filename", help="param file  (usgs format)")
parser.add_argument("dt", help="time step", type=float)
parser.add_argument("duration", help="max duration of stfs (or more)", type=float)
args = parser.parse_args()


duration = args.duration
dt = args.dt
npts = int(np.ceil(duration / dt))
print(f"duration and dt hardcoded to {duration} and {dt}s")
point_sources = {}


with open(args.filename) as fh:
    # number of segments
    line = fh.readline().strip()
    if not line.startswith("#Total number of fault_segments"):
        raise ValueError("Not a valid USGS param file.")
    nseg = int(line.split()[-1])

    # parse all segments
    for p in range(nseg):
        fault_tag = 3 if p == 0 else 64 + p
        sources = []
        axyz = []
        amt = []
        aNormMRF = []
        segment_indices = []
        # got to point source segment
        nx = 1
        for line in fh:
            # line = line.decode()
            if "#Fault_segment" in line:
                items = line.split()
                nx = int(items[4])
                ny = int(items[9])
            if "#Lat. Lon. depth" in line:
                break
        # read all point sources until reaching next segment
        for k, line in enumerate(fh):
            # line = line.decode()
            if "#Fault_segment" in line:
                break
            # Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo
            (
                lat,
                lon,
                dep,
                slip,
                rake,
                stk,
                dip,
                tinit,
                trise,
                tfall,
                M0,
            ) = map(float, line.split())

            # Negative rupture times are not supported with the current
            # logic.
            if tinit < 0:  # pragma: no cover
                raise ValueError(
                    "File contains a negative rupture time "
                    "which Instaseis cannot currently deal "
                    "with."
                )

            # Calculate the end time.
            endtime = trise + tfall

            if endtime > (npts - 1) * dt:
                raise ValueError(
                    "Rise + fall time are longer than the "
                    "total length of calculated slip. "
                    "Please use more samples."
                )

            dep *= -1e3  # km > m
            slip *= 1e-2  # cm > m
            M0 *= 1e-7  # dyn / cm > N * m
            strike = np.radians(stk)
            dip = np.radians(dip)
            rake = np.radians(rake)
            mt = compute_moment_tensor_NED(strike, dip, rake, M0)

            axyz.append([lon, lat, dep])
            amt.append(mt)
            j = k // nx
            i = k - j * nx
            segment_indices.append((i, j))
            assert npts >= (tinit + trise + tfall) / dt
            stf = asymmetric_cosine(trise, tfall, npts, dt)
            stf = np.roll(stf, round(tinit / dt))
            aNormMRF.append(stf)

    point_sources[fault_tag] = {
        "moment_tensors": cmt.NED2RTP(np.array(amt)),
        "moment_rate_functions": np.array(aNormMRF),
        "locations": np.array(axyz),
        "segment_indices": np.array(segment_indices),
    }

proj = True
potency = False
json_str = """{
    "DH": 200000.0,
    "NZ": 2,
    "slip_threshold": -10000000000.0,
    "use_geometric_center": true
}"""

bn = os.path.basename(args.filename)
cmt.write_point_source_file(f"{bn}.h5", point_sources, dt, proj, potency, json_str)
