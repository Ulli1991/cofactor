# Using cofactor

This guide goes from a first run to your own checks. Installation (clang + Enzyme) is in [INSTALL.md](INSTALL.md);
background on what to check in large simulation codes is in [HYDRO_CODES.md](HYDRO_CODES.md).

## 1. The idea in one paragraph

Numerical codes contain many hand-written mathematical claims: "this is the derivative of that", "this force belongs to
that potential", "this adjoint is the transpose of that operator", "this booking equals the energy change". cofactor
checks such claims on the **real code**: it copies the functions verbatim from your source tree, compiles them together
with a small harness under clang with the Enzyme plugin, and compares what the code claims with what automatic
differentiation of the same code gives. Every comparison is a *row* with a worst-case value and a pass threshold (*bar*)
that you write down before the first run.

The typical use: you add a physical process, write down its equations, discretise and implement them, and check with
cofactor that the implementation is that discretisation. No test problem with a known solution is needed: the extracted
functions run on chosen inputs (or on states dumped from a short real run). This complements physical tests, it does not
replace them: those tell you whether the discretisation gives the right physics; cofactor tells you whether the code is
the discretisation you meant, and points to the term and input where it is not.

## 2. First run

```
export PATH=<cofactor>/bin:$PATH
mkdir -p ~/.config/cofactor && cp <cofactor>/site.example.ini ~/.config/cofactor/site.ini   # edit it, see section 7
cofactor check-env           # finds the Enzyme plugin? which checks are visible?
cofactor list                # checks of the current project and the shipped examples
```

Run the self-contained example. On a laptop or workstation (`require_batch = no` in the site file):

```
cd <cofactor>/examples/c_cooling
cofactor run cooling
```

On a cluster where compute steps must run in jobs (`require_batch = yes`):

```
cofactor submit examples/c_cooling/cooling
```

Expected summary:

```
===== summary
  c_cooling/cooling clean           PASS
  c_cooling/cooling planted_bug     FAIL
cofactor run c_cooling/cooling: 2/2 configs as expected
```

`planted_bug` compiles a deliberately wrong derivative and is *expected* to fail (`expect = FAIL` in its config), so
"2/2 as expected" means the tool both passes correct code and catches the wrong one.

## 3. Setting up your own project

cofactor keeps a project's checks next to its code, like `.clang-format`:

```
my_code/
  src/gravity.c
  .cofactor/
    project.ini
    checks/
      softening/
        target.ini   spec.py   harness.cc   bars.txt
```

```
cd my_code
cofactor init                 # creates .cofactor/project.ini and checks/example
cofactor new softening        # adds .cofactor/checks/softening from the template
```

Inside the project, commands need no paths: the source tree defaults to the project root (the directory holding
`.cofactor/`). Commit `.cofactor/` with your code.

## 4. Writing a check, step by step

Example: `src/gravity.c` has a softened potential and the force that should be its derivative.

```c
double soft_phi(double u)   { /* potential, u = r/h */ ... }
double soft_force(double u) { /* (1/u) dphi/du, hand-coded */ ... }
```

### 4.1 `spec.py`: what to copy

```python
import sys
import extract

def spec(X, src):
    X.functions("ext_gravity.h", [
        ("src/gravity.c", r"^double soft_phi\(",   "double soft_phi(double u)"),
        ("src/gravity.c", r"^double soft_force\(", "double soft_force(double u)"),
    ])

if __name__ == "__main__":
    extract.run(spec, sys.argv)
```

- `X.functions(out, [(file, signature_regex, new_signature), ...])` copies whole functions; the body is verbatim, only
  the signature line is replaced (use it to drop attributes such as `static inline` prefixes or GPU macros, or to change
  how global arrays are passed). `{arg}` in the new signature keeps the name of a single `int` argument.
- `X.block(out, file, start_marker, stop_marker, ...)` copies lines from a unique start marker up to a stop marker:
  constants, tables, a block of a loop body. Options: `include_stop`, `start_after` (search after another unique marker),
  `expect_dropped`/`expect_closed` (preprocessor `#if/#endif` repairs you expect when a block starts or ends inside an
  `#if`), `cut_last`, `prefix`, `suffix`.
- `X.variant(src_out, new_out, [(old, new), ...])` makes a what-if copy of an extracted block with exact-once text
  replacements: a control arm with one term corrected, or instrumentation such as an iteration counter.
- `X.literal(out, text)` writes a fixed file.

Extraction stops with exit 3 if a marker is missing or not unique, a function's braces do not balance, or a preprocessor
repair differs from what you declared. Each block is hashed into a manifest; when the source changes later, the run
stops with "target changed: <block>" until you look at the change and rerun with `--accept`.

### 4.2 `harness.cc`: what to check

```cpp
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "cofactor_enzyme.h"
#include "cofactor_checks.h"
#include "ext_gravity.h"

int main(int argc, char **argv)
{
    std::string cfg = (argc > 1) ? argv[1] : "default";
    if(!getenv("COFACTOR_BARS") || !load_bars(getenv("COFACTOR_BARS"), cfg)) return 4;
    for(int i = 0; i < 4000; i++)
    {
        double u = (i + 0.5) / 4000;
        char w[64]; snprintf(w, sizeof(w), "%s u=%.6f", cfg.c_str(), u);
        double dphi = __enzyme_fwddiff<double>((void *)soft_phi, u, 1.0);   /* d phi / du by Enzyme */
        double hand = soft_force(u), ref = dphi / u;
        rec("(1) force = (1/u) dphi/du", fabs(hand - ref) / fmax(fabs(ref), 1e-300), w);
#ifndef FWD_ONLY
        double r = __enzyme_autodiff<double>((void *)soft_phi, u);          /* reverse mode, same derivative */
        rec("(2) Enzyme reverse = forward", fabs(r - dphi) / fmax(fabs(dphi), 1e-300), w);
#endif
    }
    return cofactor_report("my_code/softening");
}
```

- `__enzyme_fwddiff<double>((void *)f, x, dx, enzyme_const, c1, ...)` gives df/dx along `dx`; arguments marked
  `enzyme_const` are not differentiated. `__enzyme_autodiff<double>((void *)f, x, ...)` is reverse mode.
- `rec(name, error, where)` keeps the worst error per row name and where it happened. `relerr`/`relerr_s` give relative
  errors of vectors; NaN always counts as the worst value. `exclude(reason, where)` counts a record you could not check.
- Put reverse-mode code under `#ifndef FWD_ONLY`: if Enzyme cannot compile reverse mode for your function, cofactor
  rebuilds forward-only and the report says `REVERSE SKIPPED` (exit 2, incomplete).
- Start with a "row 0" that shows the copied code behaves like the original: bitwise equality with outputs dumped from a
  real run (section 6) is the strong form. For checks on synthetic inputs, make the copied code see what the real build
  sees: `defines_from = {src}/Config.sh` (the code's macros) and `harness_optimize` (its optimisation flags); otherwise
  a macro or flag difference between harness and build goes unnoticed.

### 4.3 `bars.txt`: pass thresholds, written before the first run

```
# <threshold> [@config-tag] <row-name prefix>
1e-12 (1) force = (1/u) dphi/du
1e-12 (2) Enzyme reverse = forward
```

The longest matching prefix wins. A bar with `@tag` applies only to configs whose name contains `tag`, and then its row
must exist. A bar whose row never appears counts as a failure (the check silently did not run); a row without a bar is
reported but not judged. Do not loosen a bar after seeing results; write down why if a row is judged by another row.

### 4.4 `target.ini`: steps and configs

```ini
[target]
description = softened gravity: force vs potential
spec = spec.py
harness = harness.cc
bars = bars.txt
harness_args = {cfg}

[config:default]
defines =

[config:double_eps]
defines = SOFTENING_VARIANT=2
```

Each `[config:NAME]` is one compile-and-run; `defines` become `-D` flags for the harness and the copied code. Add
`expect = FAIL` (or `INCOMPLETE`) for a config that must fail. Section 8 lists every key.

### 4.5 Run it

```
cofactor extract softening -c default    # text only: see what gets copied (allowed anywhere)
cofactor run softening                   # all configs; or: cofactor submit softening
```

## 5. Reading the results

Each config writes `<work>/<project>/<check>/<source-tag>/<config>/` with `gen/` (the copied code), `extract.log`,
`compile.log`, `results.txt` and the manifest. `results.txt` has one line per row:

```
(1) force = (1/u) dphi/du      n=  4000 max_rel=7.500e-01  at default u=0.951625 | FAIL (bar 1e-12)
```

`n` is how often the row was evaluated, `max_rel` the worst value, `at` where it occurred (rerun exactly that input to
investigate), then the verdict. The last line summarises PASS/FAIL counts, exclusions and whether reverse mode ran.

Exit codes: 0 pass, 1 fail, 2 incomplete (reverse mode not built), 3 extraction error, 4 setup error.
`cofactor compare OLD NEW` compares two results files row by row (useful after porting or refactoring a harness).

A failing row is a candidate, not yet a bug. First rule out the harness: does row 0 hold? Did the extraction copy the
right block? Does the check use the same inputs as the code (a missed overwrite elsewhere in the code is a classic)?
For a fix, a control arm with `X.variant` shows whether the suspected term is the whole story.

## 6. Checks on states from a real run

For routines whose inputs are hard to synthesise (a hydro face, a feedback event), dump the states from a real run:

1. Add a test-only dump hook to the code (a compile flag such as `TEST_..._DUMP_STATES` that writes the input and output
   structs at the call site). Keep it as `hook.patch` in the check, or commit the flag to the code.
2. In `target.ini`, a `build` recipe makes an isolated build copy with the hook, a `dump` recipe runs it, and
   `extract_from = build` copies the code from that build:

```ini
[target]
extract_from = build
build =
    rsync -a --exclude .git {src}/ {build}/
    cd {build} && patch -p1 < {target}/hook.patch && make -j {make_jobs}
dump =
    {mpirun} {build}/mycode {target}/run.params
dump_output = output/state_dump.bin
harness_flags = -I{build}/include
harness_args = {dumprun}/output/state_dump.bin {cfg}
```

3. Build without FMA contraction (`-ffp-contract=off`) so the harness reproduces the code bit for bit, and make row 0
   compare the harness with the dumped outputs before any derivative row counts.

`examples/gizmo/.cofactor/checks/mhd_face` is a complete example (GIZMO's MHD face, using `adapters/gizmo/build.sh`).

## 7. Site configuration (`~/.config/cofactor/site.ini`)

| section | key | meaning |
|---|---|---|
| `[toolchain]` | `prelude` | shell lines before every step (module loads, environment) |
| | `harness_prelude` | extra lines before compiling and running a harness (e.g. `module load clang/18`) |
| | `clang`, `enzyme_plugin` | compiler and `ClangEnzyme-<N>.so` |
| `[vars]` | any | machine placeholders for recipes (`make_jobs`, `mpirun`, `mpi_cflags`, ...); override project values |
| `[batch]` | `require_batch` | `yes`: compute steps only inside a batch job; use `cofactor submit` |
| | `job_env`, `submit`, `header` | batch variable that marks a job, submit command, job-script header lines |
| `[paths]` | `work` | where builds, dumps, extractions and results go |

`cofactor submit` writes a job that runs frozen copies of cofactor and of the check, so editing either while the job
waits in the queue does not change what runs.

## 8. Reference

**CLI**: `init [-C DIR]`, `new NAME`, `check-env`, `list`, `extract|build|dump|harness CHECK -c CFG [-s SRC] [--accept]`,
`run|submit CHECK [-c CFG ...] [-s SRC] [--accept]`, `compare OLD NEW`. CHECK is a check name in the current project, a
directory containing `target.ini`, a name under `$COFACTOR_TARGETS`, or a shipped example (`examples/<project>/<check>`).

**`target.ini` keys** (in `[target]`, overridable per `[config:NAME]`; defaults from the project's `[defaults]`):
`description`, `spec`, `harness`, `bars`, `std` (C++ standard of the harness, default `c++14`), `extract_from`
(`src`|`build`), `default_src`, `build`, `build_timeout`, `build_output`, `dump`, `dump_timeout`, `dump_output`,
`harness_flags`, `harness_args`, `harness_timeout`, `defines`, `defines_from` (a Config.sh-style file, e.g.
`{build}/Config.sh`: its `NAME` / `NAME=VALUE` lines become `-D` flags, so the copied code sees the build's macros),
`harness_optimize` (default `-O2`; set it to the build's optimisation flags; `-ffp-contract=off` is always added),
`expect` (`PASS`|`FAIL`|`INCOMPLETE`); `[vars]` adds placeholders.

**Placeholders** in recipes, `harness_flags`, `harness_args`: `{src}` `{build}` `{dumprun}` `{gen}` `{work}` `{target}`
`{project}` `{cfg}` `{cofactor}` `{lib}` `{plugin}` `{clang}`, every `[vars]` key (project < check < site), every key of
the config section; `{name:q}` gives the value's words shell-quoted. `${VAR}` is left to the shell; any other unknown
`{name}` is an error.

**`project.ini`**: `[project] name`, `[defaults]` (keys every check inherits), `[vars]`.

## 9. Troubleshooting

| symptom | cause and fix |
|---|---|
| `EXTRACTION ERROR: ... 0 matches` / `2 matches` | the marker moved or is not unique: pick a more specific marker |
| `target changed: <block>` | the source changed since the last run: review the change, rerun with `--accept` |
| `preprocessor repair ... spec expects ...` | the block starts or ends inside an `#if`: set `expect_dropped`/`expect_closed` |
| `full build failed ... retrying forward-only` | Enzyme could not build reverse mode; rows under `#ifndef FWD_ONLY` are skipped, exit 2 |
| harness compile error in copied code | missing types or macros: provide them in the harness before including `ext_*.h` (keep it minimal and say so in a comment), or use an adapter (`adapters/gizmo`, `adapters/arepo`) |
| `call to 'pow' is ambiguous` | another library overloads `pow`: write `std::pow(x, (double)n)` in the harness |
| `'build' compiles or runs code: only inside a batch job` | the site requires batch jobs: use `cofactor submit` |
| a row fails by exactly 1.0 or 2.0 | often a scale problem (dividing by something that is ~0, e.g. a vanishing momentum); check the magnitudes before suspecting the code |
