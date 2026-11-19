# cofactor

Checks the real GIZMO/Stella code against itself with [Enzyme](https://enzyme.mit.edu) automatic differentiation:
Jacobians, adjoints, energy and momentum bookings, Galilean invariance, table derivatives. The code under test is never
re-typed: it is copied verbatim from a build copy, replayed on states dumped by that same build, and every derivative row
counts only after the replay reproduces GIZMO.

## How a target works
1. **Hook** (`targets/<t>/hook.patch`): a test-only GIZMO flag dumps raw states at one call site (one rank, no FMA build).
2. **Extraction** (`targets/<t>/spec.py`, `lib/extract.py`): code blocks and functions copied by unique text markers.
   Fails (exit 3) on a missing or duplicated marker, an unbalanced function, an undeclared preprocessor repair, or a block
   whose hash differs from the target's manifest ("target changed"; `--accept` after review).
3. **Harness** (`targets/<t>/harness.cc`, `lib/cofactor_glue.h`, `lib/cofactor_checks.h`): row 0 = harness reproduces
   GIZMO (records that fail it are excluded from all derivative rows); then Enzyme forward/reverse checks.
4. **Bars** (`targets/<t>/bars.txt`): set before a run, `<threshold> [@config-tag] <row prefix>`. Report: every row with
   PASS/FAIL, exclusions counted, `REVERSE SKIPPED` when the reverse-mode build failed. Exit 0 pass, 1 fail, 2 incomplete.

## Targets
| target | covers | status |
|---|---|---|
| `mhd_face` | MFM/MFV face, Dedner GLM and psi/c_h cleaning, Powell terms, final per-cell bookings | port of the 2026-11-18 MHD audit; porting check against its results |
| `mhd_cg_gradient` | constrained-gradient (CG) passes of `gradients.cc` | planned (design: `docs/DESIGN.md`) |

Planned ports: MHD gradient/kernel and KVP harnesses, Stella feedback events.

## Running
Everything runs in Slurm jobs (`run/*.sbatch`); nothing is compiled or run on login nodes. Toolchain: clang 18 with
ClangEnzyme-18 (Enzyme v0.0.302), GIZMO builds with `-ffp-contract=off`.

## Known Enzyme limits (workarounds in the harnesses)
- by-value structs as shadows alias (pass element by element);
- reverse mode on a whole hydro face asserts in `GradientUtils::unwrapM` (forward columns used; reported as REVERSE SKIPPED);
- no rule for `fflush`/`printf` (registered inactive; GIZMO warning/abort paths stubbed and the record excluded).
