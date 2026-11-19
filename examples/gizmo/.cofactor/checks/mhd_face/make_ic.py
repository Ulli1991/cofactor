"""Small MHD box for the face-dump runs: 16^3 jittered periodic box (test/mhd_kernel_volume_energy generator),
divergent B blob (test/mhd_psi_cleaning_identity), plus colliding flows in x (shock at L/2, rarefaction at 0)
and a B_y reversal across x = L/2, so that the dumped faces include shocks and active limiters."""
import sys
import numpy as np
import h5py
sys.path.insert(0, sys.argv[2])
from make_mhd_box_ics import write_ics, L  # noqa: E402

fn = sys.argv[1]
write_ics(fn, 16)
with h5py.File(fn, "r+") as f:
    x = f["PartType0/Coordinates"][:]
    d = x - L / 2
    B0 = 6.61e-6
    B = np.zeros_like(x)
    B[:, 2] = B0
    B += 0.5 * B0 * d / 0.010 * np.exp(-np.sum(d ** 2, 1) / 0.010 ** 2)[:, None]
    w = L / 32
    B[:, 1] += 0.75 * B0 * np.tanh(d[:, 0] / w)
    f["PartType0"].create_dataset("MagneticField", data=B)
    v = f["PartType0/Velocities"][:]
    v[:, 0] = -40.0 * np.tanh(5 * np.sin(2 * np.pi * d[:, 0] / L))
    f["PartType0/Velocities"][...] = v
