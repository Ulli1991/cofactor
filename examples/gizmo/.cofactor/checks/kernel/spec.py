"""example_kernel: GIZMO's SPH/MFM kernel (mesh/kernel.h), copied verbatim from a GIZMO source tree (no build needed).
usage: python3 spec.py <gizmo_source> <out_dir> <manifest_dir> [--accept]"""
import sys

import extract  # cofactor lib/ on PYTHONPATH (the cofactor CLI sets it)

def spec(X, B):
    K = "mesh/kernel.h" if (B / "mesh/kernel.h").is_file() else "src/mesh/kernel.h"  # starforge_dev / gpu_kokkos_unified layout
    # the kernel definitions (shape, normalisation, neighbour numbers) for every KERNEL_FUNCTION and NUMDIMS
    X.block("ext_kernel_defs.h", K, "#if !defined(KERNEL_FUNCTION)", "#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD")
    # the two functions under test, signatures unchanged
    X.functions("ext_kernel_fun.h", [
        (K, r"^(GIZMO_GPU_FUNCTION )?static inline void kernel_hinv\(", "static inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4)"),
        (K, r"^(GIZMO_GPU_FUNCTION )?static inline void kernel_main\(", "static inline void kernel_main(double u, double hinv3, double hinv4, double *wk, double *dwk, int mode)"),
    ])


if __name__ == "__main__":
    extract.run(spec, sys.argv)
