"""GIZMO softened-gravity kernel (mesh/kernel.h kernel_gravity) and the SPH kernel of the same shape, verbatim."""
import sys

import extract


def spec(X, B):
    K = "mesh/kernel.h" if (B / "mesh/kernel.h").is_file() else "src/mesh/kernel.h"  # starforge_dev / gpu_kokkos_unified layout
    X.block("ext_kernel_defs.h", K, "#if !defined(KERNEL_FUNCTION)", "#ifdef HYDRO_PARTITION_UNITY_IMPROVE_FD")
    X.functions("ext_kernel_fun.h", [
        (K, r"^(GIZMO_GPU_FUNCTION )?static inline void kernel_hinv\(", "static inline void kernel_hinv(double h, double *hinv, double *hinv3, double *hinv4)"),
        (K, r"^(GIZMO_GPU_FUNCTION )?static inline void kernel_main\(", "static inline void kernel_main(double u, double hinv3, double hinv4, double *wk, double *dwk, int mode)"),
        (K, r"^(GIZMO_GPU_FUNCTION )?static inline double kernel_gravity\(", "static inline double kernel_gravity(double u, double hinv, double hinv3, int mode)"),
    ])


if __name__ == "__main__":
    extract.run(spec, sys.argv)
