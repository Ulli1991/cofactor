"""c_cooling: what to copy from the code under test (here cooling.c at the project root), verbatim, by text markers.
usage (the cofactor CLI calls it): PYTHONPATH=<cofactor>/lib python3 spec.py <src> <out_dir> <manifest_dir> [--accept]"""
import sys

import extract

C = "cooling.c"


def spec(X, src):
    # the table (file-scope constants the functions read)
    X.block("ext_table.h", C, "static const double T_edges[5]", "/* Lambda(T) [erg cm^3 s^-1]")
    # the functions under test, signatures unchanged
    X.functions("ext_cooling.h", [
        (C, r"^double cooling_lambda\(", "double cooling_lambda(double T)"),
        (C, r"^double cooling_dlambda_dT\(", "double cooling_dlambda_dT(double T)"),
        (C, r"^double cooling_residual\(", "double cooling_residual(double T, double T0, double n2dt)"),
        (C, r"^double cooling_residual_dT\(", "double cooling_residual_dT(double T, double T0, double n2dt)"),
    ])


if __name__ == "__main__":
    extract.run(spec, sys.argv)
