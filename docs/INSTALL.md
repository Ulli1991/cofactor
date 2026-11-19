# Installing cofactor

cofactor itself is Python 3.6+ (standard library only) and a few C++ headers; nothing to compile. It needs a clang with the
Enzyme plugin for the same LLVM major version.

## 1. clang + Enzyme
Tested: clang/LLVM 18.1.8, Enzyme commit `b1d93d0b` (v0.0.302 era), Release build, plugin `ClangEnzyme-18.so`.

    git clone https://github.com/EnzymeAD/Enzyme && cd Enzyme/enzyme
    mkdir build && cd build
    cmake .. -DLLVM_DIR=<llvm-18>/lib/cmake/llvm -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=<llvm-18>/bin/clang -DCMAKE_CXX_COMPILER=<llvm-18>/bin/clang++
    make -j 8 ClangEnzyme-18          # -> build/Enzyme/ClangEnzyme-18.so

Build it with the same clang you will compile harnesses with (on a cluster: inside a batch job). If your system has
clang as a module, `-DLLVM_DIR=$(dirname $(which clang))/../lib/cmake/llvm` usually works.

## 2. cofactor
    git clone https://github.com/Ulli1991/cofactor
    export PATH=$PWD/cofactor/bin:$PATH
    mkdir -p ~/.config/cofactor && cp cofactor/site.example.ini ~/.config/cofactor/site.ini   # then edit it
    cofactor check-env

`site.ini` describes the machine: module loads (`prelude`), the clang and plugin paths, placeholders such as `mpirun`,
whether compute steps must run in a batch job (`require_batch`) and the batch header for `cofactor submit`.

## 3. Smoke test
    cofactor submit examples/c_cooling/cooling      # or: cd cofactor/examples/c_cooling && cofactor run cooling
Expected: `clean` PASS, `planted_bug` FAIL (a deliberately wrong derivative must be caught).

## Notes
- Harnesses are compiled with `-ffp-contract=off`; build the code under test the same way when the harness has to
  reproduce its results bit for bit.
- If the full (reverse-mode) compile fails inside Enzyme, cofactor retries forward-only and reports
  `REVERSE SKIPPED` (exit 2, incomplete) instead of hiding it.
