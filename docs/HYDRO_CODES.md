# cofactor with large hydrodynamics codes

Large astrophysical codes (GIZMO, AREPO, GADGET, Athena++, ENZO, FLASH, RAMSES, ...) contain many hand-written derivatives,
linearisations and bookkeeping identities that are rarely tested directly: Riemann-solver states, slope limiters, gradient
estimators, implicit cooling and chemistry Jacobians, divergence-cleaning energy exchanges, feedback deposition weights.
They are usually validated only through end-to-end tests (Sod, Orszag-Tang, a cloud collapse), where a small error in one
term is hidden by everything else. cofactor checks such pieces one at a time, in the real code, with automatic
differentiation as the reference.

## What it is for

| question | how cofactor answers it | example |
|---|---|---|
| Is a hand-coded derivative right? | Enzyme forward mode of the same function | kernel `dwk`, cooling `dLambda/dT`, Newton Jacobians of implicit solvers |
| Is a hand-coded adjoint/transpose right? | Enzyme Jacobian vs the code's adjoint; identity <Ju, v> = <u, J^T v> | psi/c_h cleaning adjoint of the face divergence (GIZMO), gradient-estimator transpose |
| Does an energy (or mass, momentum) booking close? | Enzyme derivative of the real energy function along the code's own rates = the booked rate | MHD energy ledger rows, Powell terms, upwind heat, feedback energy/momentum per event |
| Is a scheme invariant under a symmetry? | derivative with respect to the symmetry parameter = 0 | Galilean boost of a face flux (GIZMO MFM/MFV; Powell-into-heat variants were shown to break it exactly by 2 D B_k) |
| Is a function homogeneous as assumed? | Euler identity x . df/dx = f | face divergence in the current B (fails when gradients carry history) |
| Is a table interpolation smooth where it must be? | derivative vs central differences inside cells, jump scans at cell edges | stellar yield tables, cooling tables |
| Does the extracted/re-implemented code still equal the original? | row 0: harness = code, bitwise on dumped states | every check |

## How it works on a big code

1. **Find the routine and its inputs.** Pick a function (or a block of a loop body) with a clear question, e.g. "the
   per-face flux and the psi adjoint it books". Note every input it reads: function arguments, the particle/cell arrays,
   global parameters, neighbour data.
2. **Dump real states (hook).** Add a test-only compile flag in the code (e.g. `TEST_..._DUMP_STATES`) that writes the raw
   input and output structs at that call site, on one MPI rank, for a sample of calls. Keep it as `hook.patch` in the check
   (or commit the flag to the code). Real states contain the limiter branches, fallbacks and edge cases synthetic inputs miss.
3. **Build without FMA contraction** (`-ffp-contract=off`) and with the hook on. The harness, compiled by clang, can then
   reproduce the code's arithmetic exactly; row 0 is bitwise, so any later difference is a real derivative difference.
4. **Copy the code, don't re-type it.** `spec.py` copies the routine verbatim by text markers from the build copy.
   Functions that read global particle arrays get a new signature taking the arrays as arguments; inside the harness the
   global names are macros for the arguments (Enzyme cannot shadow pointer globals). Abort and warning paths are stubbed
   and counted (`adapters/gizmo/gizmo_glue.h` is the pattern); records that reach them are excluded and counted.
5. **Replay and differentiate.** The harness reads each dumped record, re-runs the copied code (row 0), then asks Enzyme
   for the derivatives the question needs: columns of a Jacobian in forward mode, adjoints in reverse mode, derivatives
   along the code's own rates for bookings, along a symmetry for invariance.
6. **Bars, report, record.** Bars are fixed before the run. The report lists every row, PASS/FAIL, the worst record
   (file position, step, face/particle) so a failure can be replayed alone, exclusions and whether reverse mode ran.

## Practical points (from the GIZMO audits)
- One rank, OpenMP off, a small box (16^3) that still reaches shocks, limiters and fallbacks; a few thousand records
  suffice. Sample from several steps; the first step is often special.
- Dump both the state before and after the call: running sums (accelerations, dB/dt) can only be compared as
  after - before, to the rounding of the accumulated value.
- Compare per term, not only totals: a booking that closes in total can still hide two compensating errors.
- When a check fails, first rule out the harness (row 0, the extraction, a missed input such as a density overwrite
  elsewhere in the code). One GIZMO "failure" (KVP force exact only for equal masses) was the harness reading the wrong
  density; the fix was in the harness, not the code.
- Reverse mode on a whole hydro face can hit Enzyme limitations; forward mode, one input direction at a time, gives the
  same Jacobian and is what the GIZMO checks use. cofactor reports this as REVERSE SKIPPED.
- Keep the code's own checks too: cofactor complements conservation tests, convergence tests and comparisons with
  reference codes; it does not replace them.

## Codes
- **GIZMO**: adapter `adapters/gizmo` (stubs for `PRINT_WARNING`/`endrun`/`fflush`, build script writing `Config.sh`,
  no-FMA `OPTIMIZE`). Example setup: `examples/gizmo/.cofactor` (kernel derivative; MHD face with dump hook).
- **AREPO, GADGET-4**: the same pattern; AREPO's `terminate`/`warn` macros are the abort/warning paths to stub, `Config.sh`
  and `OPT` in the Makefile take the role of GIZMO's `Config.sh`/`OPTIMIZE`. Keep the setup in the code's own repository
  (`.cofactor/`), especially for codes that are not public: a check's spec, harness and hook contain pieces of the code.
- **Fortran codes** (RAMSES, FLASH, ...): Enzyme works on LLVM IR, so Fortran via flang is possible in principle; cofactor's
  harnesses are C++, calling Fortran routines through `bind(C)` interfaces. Not tested yet.
- **GPU codes** (Kokkos, CUDA): check the CPU path with cofactor and compare the GPU kernels against it numerically.
