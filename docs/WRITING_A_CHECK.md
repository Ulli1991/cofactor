# Writing a check

A check answers one question about real code with Enzyme: does a hand-coded derivative, Jacobian, adjoint or energy
booking equal what automatic differentiation of the same code gives? The code under test is never re-typed.

## Layout (in your project)
    .cofactor/project.ini              shared defaults and placeholders
    .cofactor/checks/<name>/target.ini steps and configs
    .cofactor/checks/<name>/spec.py    what to copy from the source, by text markers
    .cofactor/checks/<name>/harness.cc calls the copied code, records rows
    .cofactor/checks/<name>/bars.txt   pass thresholds, written before the first run
Start with `cofactor init` (once per project) and `cofactor new <name>`.

## The four rules
1. **Copy, don't re-type.** `spec.py` copies functions and blocks verbatim. Markers must be unique; a moved or edited
   block stops the run ("target changed", exit 3) until you review it and pass `--accept`.
2. **Row 0 first.** The first rows show that the harness reproduces the original (same outputs from the same inputs,
   ideally bitwise from a dump of a real run). Records that fail row 0 get no derivative rows.
3. **Bars before results.** Thresholds go into `bars.txt` before the first run and are not loosened afterwards.
   A bar whose row never appears counts as a failure.
4. **Say what was not checked.** Exclusions are counted; a forward-only build reports `REVERSE SKIPPED`.

## Three kinds of checks
- *Synthetic inputs* (examples/c_cooling, examples/gizmo kernel): copy functions, call them on chosen inputs.
- *Dumped states* (examples/gizmo mhd_face): a test-only hook in the code (`hook.patch`, applied by the build recipe)
  writes raw states from a real run; the harness replays them. Needs `build`, `dump`, `extract_from = build`.
- *Identities*: Euler/homogeneity (x . df/dx = f), adjoint pairs (<Ju, v> = <u, J^T v>), conservation of a booked
  quantity along the code's own rates, invariance under a symmetry (e.g. a Galilean boost: d/dV = 0).

## Rows
`rec("name", error, "where")` keeps the worst error per row name and where it happened. Use `relerr`/`relerr_s` for
relative errors; NaN always fails. Name rows so that bar prefixes select them; tag bars with `@text` to require them
only for configs whose name contains `text`.

## Enzyme tips
- Pass structs element by element or as arrays; by-value struct shadows can alias.
- Globals that are not differentiated: `__attribute__((used)) void *__enzyme_inactive_global_X = (void *)&X;`
- Library calls without a derivative rule (printf, fflush, abort paths): register inactive or stub them (see
  adapters/gizmo/gizmo_glue.h).
- If reverse mode on a large function fails, forward mode column by column gives the same Jacobian.
