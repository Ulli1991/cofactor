# cofactor

**cofactor is a framework for checking numerical code against itself with automatic differentiation.**
It takes the real code under test (copied verbatim, never re-typed), runs it on real or synthetic inputs, and uses
[Enzyme](https://enzyme.mit.edu) to verify what the code claims about itself: hand-coded derivatives and Jacobians,
adjoints and transposes, conservation bookings, invariances. It was built for large simulation codes, where such pieces
are numerous and usually tested only end to end, but it works on any C/C++ code that clang compiles.

Like `.clang-format` or `.clang-tidy`, a project keeps its own small setup in a `.cofactor/` directory; cofactor is the tool
that runs it.

```
my_code/
  src/...
  .cofactor/
    project.ini                 shared defaults and placeholders
    checks/
      flux_jacobian/
        target.ini              steps (build, dump, extract, harness) and configs
        spec.py                 what to copy from the source, by unique text markers
        harness.cc              calls the copied code, records check rows
        bars.txt                pass thresholds, fixed before the first run
```

## Quick start
```
export PATH=<cofactor>/bin:$PATH
cp <cofactor>/site.example.ini ~/.config/cofactor/site.ini   # machine: clang + Enzyme plugin, modules, batch system
cofactor check-env
cd <cofactor>/examples/c_cooling && cofactor run cooling     # or 'cofactor submit cooling' on a cluster
```
Expected: `clean` PASS, `planted_bug` FAIL (the example plants a wrong derivative to show it is caught).

For your code: `cofactor init` at the project root, `cofactor new <check>`, edit the four files, `cofactor run <check>`.

## What a run does
1. **build** (optional recipe): an isolated build copy, e.g. with a test-only hook that dumps raw states from a real run;
2. **dump** (optional recipe): the run that writes those states;
3. **extract**: `spec.py` copies functions/blocks verbatim. A missing or duplicated marker, an unbalanced function, or a
   block that changed since the last run stops here (exit 3) until reviewed (`--accept`);
4. **harness**: compiled with clang + Enzyme (`-ffp-contract=off`); row 0 shows the harness reproduces the code, then the
   derivative/identity rows; if Enzyme's reverse mode fails to compile, the forward-only build runs and the report says
   `REVERSE SKIPPED`;
5. **report**: every row with its worst value and where it occurred, PASS/FAIL against `bars.txt`, exclusions counted.
   Exit 0 pass, 1 fail, 2 incomplete, 3 extraction error, 4 setup error.

On clusters, compute steps can be restricted to batch jobs (`require_batch` in the site file); `cofactor submit` writes and
submits the job, which runs frozen copies of cofactor and of the check.

## Documentation
- [docs/USAGE.md](docs/USAGE.md): **how to use it**: first run, setting up a project, writing a check step by step,
  reading results, checks on states from real runs, configuration reference, troubleshooting
- [docs/INSTALL.md](docs/INSTALL.md): clang + Enzyme, site configuration, smoke test
- [docs/WRITING_A_CHECK.md](docs/WRITING_A_CHECK.md): the four files, the rules, kinds of checks, Enzyme tips
- [docs/HYDRO_CODES.md](docs/HYDRO_CODES.md): using cofactor with large hydro codes (GIZMO, AREPO, ...) and what for

## Examples
| project | check | shows |
|---|---|---|
| `examples/c_cooling` | `cooling` | the smallest complete check: plain C, synthetic inputs, a planted bug that must fail |
| `examples/gizmo` | `kernel` | real code (GIZMO's SPH kernel), every kernel shape x dimension, no build needed |
| `examples/gizmo` | `mhd_face` | a large code end to end: hook, no-FMA build, dump run, bitwise replay of the MHD face, adjoints, energy bookings, Galilean invariance |

The GIZMO setup is run against a GIZMO source tree: `cofactor run examples/gizmo/kernel -s <gizmo>`.

## Status
Early. Tested with clang/LLVM 18.1.8 and Enzyme `b1d93d0b` on Linux with Slurm. C/C++ code; Fortran via flang untested.
