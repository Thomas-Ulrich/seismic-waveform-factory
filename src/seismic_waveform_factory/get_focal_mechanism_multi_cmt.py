#!/usr/bin/env python3
import argparse

import h5py
from obspy.imaging.beachball import MomentTensor, mt2plane

parser = argparse.ArgumentParser(
    description="compute nodal plate from Moment Tensor (h5) file"
)
parser.add_argument("filename", help="point source file (h5)")
args = parser.parse_args()

h5f = h5py.File(args.filename, "r")
aMomentTensor = h5f["moment_tensors"][:, :]
nsource = aMomentTensor.shape[0]
print("source\tstrike\tdip\trake")
for i in range(nsource):
    mt = MomentTensor(aMomentTensor[i, :], 1)
    nodalplane = mt2plane(mt)
    print(f"{i}\t{nodalplane.strike}\t{nodalplane.dip}\t{nodalplane.rake}")
