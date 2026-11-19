# Design: Enzyme debugger package for GIZMO/Stella - draft 2026-11-19

Goal (Uli 2026-11-19): generalise the two Enzyme audits (MHD faces/gradients/KVP; Stella feedback events) into one reusable
tool that checks any GIZMO routine's derivatives, Jacobians, adjoints and conservation bookings against the real code, bit for
bit; later possibly a CS paper. First new customer: the psi/c_h constrained-gradient (CG) measurement.

## What exists (to be generalised, not rewritten)
| piece | MHD audit (enzyme/mhd_audit) | feedback audit (enzyme/fb_audit) |
|---|---|---|
| dump hook in GIZMO (test-only flag) | TEST_MHD_DUMP_FACE_STATES, hydro_evaluate.h (commit 1270c3c7b, patch file) | TEST_FB_DUMP_EVENT_STATES (patch_dump.py, commit c6d6d68a2) |
| extractor: verbatim code by text markers | harness/gen.py, gen2.py | harness/gen.py (struct_block, function with new signature) |
| harness: replay = GIZMO, Enzyme fwd/rev, check rows | face_harness.cc, grad_harness.cc, kvp_harness.cc | fb_harness.cc |
| build/run | job*_*.sh, harness/build_and_run.sh (clang/18 + ClangEnzyme-18.so, FWD_ONLY fallback) | job1/2/3.sh, harness/build_and_run.sh |
Common conventions already in use: no-FMA builds (-ffp-contract=off) so the harness reproduces GIZMO exactly; P/CellP passed
as arrays (macros); stubs for fflush/PRINT_WARNING/endrun; inactive globals (All); check rows "name n= max_rel= at <where>";
row 0 "harness = GIZMO run" must be bitwise/round-off before any derivative row counts.

## Package layout (proposal)
gizmo-enzyme-debug/  (own local git repo, like stella-tables; GIZMO side only gets TEST_ dump hooks behind flags)
- core/enzyme_glue.h     Enzyme declarations, inactive-function/global registrations, GIZMO stubs, P/CellP array macros
- core/checks.h          rec()/relerr()/check table, worst-case location, MISMATCH lines, branch categories, results.txt writer
- core/dump_io.h         one record format for every hook: header (git hash, Config.sh hash, struct sizes, rank/step/time),
                         raw structs; reader validates the fingerprint against the build it is compiled with
- extract/extract.py     marker-based block/function/struct extraction with signature rewriting; per-target spec file
- targets/<name>/        spec.py (what to extract), hook.patch (GIZMO dump hook), harness.cc (checks), README (bars)
                         first targets: mhd_face, mhd_grad, mhd_kvp, fb_event (ported from the audits, results must match)
                         new: mhd_cg_gradient
- run/                   one sbatch template: isolated build copies (no FMA) from a worktree + hook patch, dump runs with
                         process-group timeout, harness compile with fwd-only fallback, collation; nothing on login nodes
- docs/LIMITS.md         Enzyme limits met so far (by-value struct shadows alias, reverse mode on full face asserts in
                         GradientUtils::unwrapM, no fflush rule) and the workarounds

## First new target: mhd_cg_gradient (answers the open CG questions)
Hook: dump the full cell state of every active cell and its neighbours at entry to the gradient routine and after each CG pass
(one rank, Dedner+CG default build, box and mixed-timebin runs). Harness: replay the CG passes (gradients.cc pair code and the
between-pass block, verbatim) and measure, per face:
 (1) history term h = D(current B, stale G) - D(current B, G from a fresh start), relative to |D| and to psi V D;
 (2) inactive-neighbour remainder (mixed timebins);
 (3) Euler residual B.dD/dB - D with Enzyme forward mode through all passes (frozen and unfrozen limiter factors);
 (4) Enzyme Jacobian of the full CG map, the reference for a later hand adjoint (option B).
Bars to be fixed before the run; result decides option A/B in docs/DESIGN_psi_constrained_gradient.md.

## Porting check
Each ported target must reproduce its audit's results.txt (same rows, same max_rel) on the same dumps before the old
harness is retired.

## Open questions for Uli
- Repo location/name (proposal: /u/uli/gizmo-enzyme-debug, local git, not pushed).
- Kokkos/GPU path: harness on the CPU code only; GPU kernels checked against the Enzyme-verified CPU path (as decided).

## Review (fresh-context agent, 2026-11-19) and revised plan
Findings adopted (spot-checked: FWD_ONLY fallback leaves no row, build_and_run.sh:16-23; extraction stops at the dump-hook
marker TEST_FACE_DUMP_B, gen.py:73):
- Fingerprint today is sizeof only (both dumps). New: magic/version, struct sizes + offsetof of key fields, git hash,
  Config.sh hash, compile flags, hash of every extracted block; a record whose fingerprint differs is rejected, not printed.
- Extraction fails loudly: missing or duplicate marker, brace-unbalanced span (brace matching, not "first line starting
  with }"), any balance() repair -> error. EXTRACTED.txt stores block hashes; a changed block is reported as "target changed".
- Loud fallback: a REVERSE SKIPPED row whenever the fwd-only build is used. Exclusions (abort/warn paths, states not
  reproduced) counted in every report. Row 0 (harness = GIZMO) gates: derivative rows are not evaluated if it fails.
- Bars live in each target's spec and are printed PASS/FAIL; set before the run.
- Hook patches: a check that each patch applies and the dump size/fingerprint is unchanged, run before every audit.
Scope cut (over-engineering): no generic spec language, no single payload format (shared header + per-target payload),
no LIMITS/paper scaffolding yet.

Minimum viable first version:
1. lib: checks.h (rows, bars, exclusions, REVERSE SKIPPED, row-0 gate), fingerprint.h, extract.py (block/function with the
   failure rules above). Port mhd_face only; it must reproduce results.txt of 231637/231830 row by row on the same dumps.
2. mhd_cg_gradient, forward mode only (reverse on a whole pass is not realistic, cf. GradientUtils assert):
   - routine: hydro/gradients.cc, the gradient_iteration loop (674) with the per-pair CG block (1560-1642) and the
     between-pass block (887-1029). Hook TEST_MHD_DUMP_CG_GRADIENT_STATES: at entry and after each pass, for sampled active
     cells: P/CellP of the cell and of every neighbour (pairs as in the gradient loop), Gradients.B before the pass,
     GasGradDataPasser fields; one rank. Size: ~2 kB per cell-neighbour record; cap 200 cells/step x 10 steps.
   - builds: Dedner + MHD_CONSTRAINED_GRADIENT (=1), MFM and MFV, no FMA; runs: the face-dump box (one timebin) and the
     mixed-timebin identity setup.
   - quantities and bars (fixed now):
     row 0: replay of the CG passes = GIZMO final Gradients.B, relative <= 1e-14 (bitwise expected);
     (1) history term per face h_f = D_f(B, G_stale) - D_f(B, G_fresh), G_fresh = same replay with the pass-0 neighbour
         gradients replaced by this step's unconstrained ones (that replay deliberately differs from the run); reported as
         |h_f|/|D_f| (median, 99 %, max) and sum_f |psi_f h_f| / sum_f |psi_f D_f| with psi from the dump; informational,
         decides A vs B (A requires the latter <= 1e-12, i.e. negligible; otherwise B);
     (2) inactive-neighbour remainder: same ratio restricted to faces whose neighbour was inactive on the step (mixed runs);
     (3) Euler residual per face |B.dD_f/dB - D_f| / (|B_i| A_f): with frozen factors (expected = |h_f|-part only) and with
         Enzyme through all factors; bar for option B: <= 1e-12 with fresh start;
     (4) forward-mode Jacobian columns of the full CG map, stored as the reference for the hand adjoint (no bar yet).
3. Repo: /u/uli/gizmo-enzyme-debug (local git, not pushed) pending Uli's answer.
