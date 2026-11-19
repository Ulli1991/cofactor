"""What to copy from the code under test, verbatim, by text markers (cofactor lib/extract.py).
Run by the cofactor CLI: PYTHONPATH=<cofactor>/lib python3 spec.py <src> <out_dir> <manifest_dir> [--accept]"""
import sys

import extract


def spec(X, src):
    # X.block(out_name, file, start_marker, stop_marker, include_stop=False, start_after=None,
    #         expect_dropped=0, expect_closed=0)            -> lines from the unique start marker to the stop marker
    # X.functions(out_name, [(file, signature_regex, new_signature), ...])  -> whole functions, bodies verbatim
    X.functions("ext_code.h", [
        ("TODO/path/to/file.c", r"^double my_function\(", "double my_function(double x)"),
    ])


if __name__ == "__main__":
    extract.run(spec, sys.argv)
