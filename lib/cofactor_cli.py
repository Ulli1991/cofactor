"""cofactor: a framework for checking numerical code against itself with Enzyme automatic differentiation.

A project keeps its checks in a .cofactor/ directory at its root (like .clang-format / .clang-tidy):
    <project>/.cofactor/project.ini          project-wide defaults ([defaults]) and placeholders ([vars])
    <project>/.cofactor/checks/<name>/       one check: target.ini, spec.py, harness.cc, bars.txt
cofactor finds the project by walking up from the current directory (or -C DIR).

  cofactor init [-C DIR]                              create .cofactor/ in a project (with one template check)
  cofactor new NAME                                   add a check to the current project from the template
  cofactor check-env                                  site configuration and toolchain
  cofactor list                                       checks of the current project and the shipped examples
  cofactor extract CHECK [-s SRC] -c CFG [--accept]   copy the code under test (text only; allowed anywhere)
  cofactor build   CHECK [-s SRC] -c CFG              the check's build recipe (batch only)
  cofactor dump    CHECK [-s SRC] -c CFG              the check's dump recipe (batch only)
  cofactor harness CHECK [-s SRC] -c CFG              extract, compile and run the harness, report (batch only)
  cofactor run     CHECK [-s SRC] [-c CFG ...]        all steps for each config, summary (batch only)
  cofactor submit  CHECK [-s SRC] [-c CFG ...]        batch job running 'cofactor run' on frozen copies of cofactor and the check
  cofactor compare OLD NEW                            two results files, row by row

CHECK: a name in the current project, a directory with target.ini, a name in $COFACTOR_TARGETS (':'-separated), or a
shipped example (examples/<project>/<check>). SRC: the source tree under test; default: the project root (or the check's
default_src). Site configuration (machine, toolchain, batch system): $COFACTOR_SITE, ~/.config/cofactor/site.ini, or
<cofactor>/site.ini (template: site.example.ini). Exit codes: 0 pass, 1 fail, 2 incomplete (reverse mode not built),
3 extraction error, 4 setup error.
"""
import argparse
import configparser
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORD = {0: "PASS", 1: "FAIL", 2: "INCOMPLETE", 3: "EXTRACTION ERROR", 4: "SETUP ERROR"}


class Fail(Exception):
    def __init__(self, msg, code=4):
        super().__init__(msg)
        self.code = code


def ini(path):
    c = configparser.ConfigParser(interpolation=None)
    c.optionxform = str
    c.read(path)
    return c


def lines(v):
    return [l for l in (v or "").split("\n") if l.strip()]


# ---------------------------------------------------------------- site, project, checks
START = Path.cwd()


def find_project(start):
    """the nearest directory at or above start that holds .cofactor/"""
    p = Path(start).resolve()
    for d in [p] + list(p.parents):
        if (d / ".cofactor").is_dir():
            return d
    return None


def load_site():
    cands = [os.environ.get("COFACTOR_SITE"), Path.home() / ".config" / "cofactor" / "site.ini", ROOT / "site.ini"]
    for c in cands:
        if c and Path(c).is_file():
            S = ini(c)
            S.path = Path(c)
            return S
    raise Fail("no site configuration: copy site.example.ini to ~/.config/cofactor/site.ini (or <cofactor>/site.ini, "
               "or set COFACTOR_SITE) and edit it")


def check_sources():
    """(label, directory holding checks) in search order"""
    out = []
    P = find_project(START)
    if P is not None:
        out.append((P.name, P / ".cofactor" / "checks"))
    for x in os.environ.get("COFACTOR_TARGETS", "").split(":"):
        if x:
            out.append((x, Path(x)))
    ex = ROOT / "examples"
    if ex.is_dir():
        for proj in sorted(ex.iterdir()):
            if (proj / ".cofactor" / "checks").is_dir():
                out.append((f"examples/{proj.name}", proj / ".cofactor" / "checks"))
    return out


def all_targets():
    out = {}
    for label, base in check_sources():
        if base.is_dir():
            for t in sorted(base.rglob("target.ini")):
                rel = str(t.parent.relative_to(base))
                key = rel if label == (find_project(START).name if find_project(START) else None) else f"{label}/{rel}"
                out.setdefault(key, t.parent)
    return out


def load_target(name):
    d = Path(name)
    if not (d / "target.ini").is_file():
        d = all_targets().get(name)
        if d is None:
            raise Fail(f"unknown check {name!r}; 'cofactor list' shows the known ones")
    d = d.resolve()
    T = ini(d / "target.ini")
    T.dir = d
    T.name = d.name
    T.project = find_project(d)  # the project the check belongs to (None for a free-standing check directory)
    T.P = ini(T.project / ".cofactor" / "project.ini") if T.project and (T.project / ".cofactor" / "project.ini").is_file() else ini("/dev/null")
    if T.project is not None:
        T.name = f"{T.project.name}/{d.name}"
    return T


def configs(T):
    return [s.split(":", 1)[1] for s in T.sections() if s.startswith("config:")]


# ---------------------------------------------------------------- placeholders
PH = re.compile(r"(\$?)\{(\w+)(:q)?\}")


def context(S, T, cfg, src):
    W = workdir(S, T, cfg, src)
    ctx = {"cofactor": str(ROOT), "lib": str(ROOT / "lib"), "target": str(T.dir), "cfg": cfg, "work": str(W),
           "build": str(W / "build"), "dumprun": str(W / "dumprun"), "gen": str(W / "gen"), "src": str(src) if src else "",
           "plugin": os.path.expandvars(S.get("toolchain", "enzyme_plugin", fallback="")),
           "clang": S.get("toolchain", "clang", fallback="clang++")}
    ctx["project"] = str(T.project) if T.project else ""
    for sec in (T.P, T, S):  # project, then check, then this machine's site values; config keys below win over all
        if sec.has_section("vars"):
            ctx.update(dict(sec.items("vars")))
    if T.has_section(f"config:{cfg}"):
        ctx.update(dict(T.items(f"config:{cfg}")))
    return ctx


def subst(text, ctx, what="recipe"):
    """{key} -> value, {key:q} -> its words shell-quoted; ${...} left to the shell; any other {key} is an error"""
    for _ in range(4):  # values may contain placeholders
        def rep(m):
            if m.group(1):
                return m.group(0)
            k = m.group(2)
            if k not in ctx:
                raise Fail(f"{what}: unknown placeholder {{{k}}} (known: {', '.join(sorted(ctx))})")
            v = ctx[k]
            return " ".join(shlex.quote(w) for w in v.split()) if m.group(3) else v
        new = PH.sub(rep, text)
        if new == text:
            return new
        text = new
    return text


def tval(T, cfg, key, default=""):
    """config section, then the check's [target], then the project's [defaults]"""
    s = f"config:{cfg}"
    if T.has_option(s, key):
        return T.get(s, key)
    if T.has_option("target", key):
        return T.get("target", key)
    return T.P.get("defaults", key, fallback=default)


def workdir(S, T, cfg, src=None):
    """<work>/<check>/<source tree tag>/<config>: runs of the same check on different source trees never share files"""
    base = Path(os.path.expandvars(S.get("paths", "work"))) / T.name
    if src:
        import hashlib
        p = Path(src).resolve()
        base = base / f"{p.name}_{hashlib.sha1(str(p).encode()).hexdigest()[:6]}"
    return base / cfg


# ---------------------------------------------------------------- shell with the site prelude and a process-group timeout
def sh(S, script, log, timeout, cwd=None, harness=False):
    pre = lines(S.get("toolchain", "prelude", fallback=""))
    if harness:
        pre += lines(S.get("toolchain", "harness_prelude", fallback=""))
    body = "\n".join(pre) + "\nset -e\n" + script + "\n"  # site prelude (module systems) runs without -e/-u
    with open(log, "w") as f:
        p = subprocess.Popen(["bash", "-c", body], stdout=f, stderr=subprocess.STDOUT, cwd=cwd, start_new_session=True)
        try:
            return p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
            f.write(f"\ncofactor: TIMEOUT after {timeout} s, process group killed\n")
            return 124


def require_batch(S, what):
    if S.get("batch", "require_batch", fallback="no").lower() in ("yes", "true", "1"):
        if not os.environ.get(S.get("batch", "job_env", fallback="SLURM_JOB_ID")):
            raise Fail(f"'{what}' compiles or runs code: only inside a batch job on this site (require_batch = yes); use 'cofactor submit'")


def tail(path, n=10):
    try:
        return "\n".join(Path(path).read_text().split("\n")[-n:])
    except OSError:
        return ""


def resolve_src(T, cfg, src, need):
    if src:
        return Path(src).resolve()
    d = tval(T, cfg, "default_src")
    if d:
        return Path(d.replace("{target}", str(T.dir)).replace("{project}", str(T.project or ""))).resolve()
    if T.project is not None:
        return T.project
    if need:
        raise Fail(f"check {T.name} is not inside a project: give the source tree with -s SRC")
    return None


# ---------------------------------------------------------------- steps
def do_extract(S, T, cfg, src, accept=False):
    W = workdir(S, T, cfg, src)
    W.mkdir(parents=True, exist_ok=True)
    frm = tval(T, cfg, "extract_from", "src")
    base = W / "build" if frm == "build" else src
    if frm == "build" and not (W / "build").is_dir():
        raise Fail(f"{T.name} extracts from its build: run 'cofactor build' first")
    # the manifest (hashes of the extracted blocks) lives in the work directory: it tracks this site's source between runs
    cmd = [sys.executable, str(T.dir / tval(T, cfg, "spec", "spec.py")), str(base), str(W / "gen"), str(W)]
    env = dict(os.environ, PYTHONPATH=str(ROOT / "lib") + os.pathsep + os.environ.get("PYTHONPATH", ""))
    r = subprocess.run(cmd + (["--accept"] if accept else []), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, timeout=300, env=env)
    (W / "extract.log").write_text(r.stdout)
    print(r.stdout.rstrip())
    if r.returncode:
        raise Fail(f"extraction failed ({W / 'extract.log'})", 3)


def do_recipe(S, T, cfg, src, key, subdir, what):
    require_batch(S, what)
    recipe = tval(T, cfg, key)
    if not recipe.strip():
        print(f"{T.name} {cfg}: no {what} step")
        return
    W = workdir(S, T, cfg, src)
    d = W / subdir
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    script = subst(recipe, context(S, T, cfg, src), f"{T.name} {key}")
    t0 = time.time()
    rc = sh(S, script, W / f"{what}.log", timeout=int(tval(T, cfg, f"{key}_timeout", "3600")), cwd=d)
    print(f"{what} {cfg}: exit {rc} ({time.time() - t0:.0f} s), log {W / (what + '.log')}")
    if rc:
        print(tail(W / f"{what}.log", 12))
        raise Fail(f"{what} of {cfg} failed", 4)
    out = tval(T, cfg, f"{key}_output")
    if out and not (d / subst(out, context(S, T, cfg, src))).exists():
        raise Fail(f"{what} of {cfg}: expected output {out} missing", 4)


def config_defines(path):
    """macro definitions from a Config.sh-style file: one NAME or NAME=VALUE per line, '#' and '%' start comments"""
    p = Path(path)
    if not p.is_file():
        raise Fail(f"defines_from: {p} not found", 4)
    out = []
    for l in p.read_text().split("\n"):
        l = re.split(r"[#%]", l, 1)[0].strip()
        if l:
            out.append(l.split()[0])
    return out


def do_harness(S, T, cfg, src, accept=False):
    require_batch(S, "harness")
    W = workdir(S, T, cfg, src)
    do_extract(S, T, cfg, src, accept)
    ctx = context(S, T, cfg, src)
    dl = [d for l in lines(tval(T, cfg, "defines")) for d in l.split()]
    df = tval(T, cfg, "defines_from")
    if df.strip():  # the code's own config file (e.g. GIZMO/AREPO/GADGET Config.sh): the copied code sees the build's macros
        dl = config_defines(subst(df, ctx, f"{T.name} defines_from")) + dl
    defs = " ".join(shlex.quote(f"-D{d}") for d in dl)
    flags = subst(tval(T, cfg, "harness_flags"), ctx, f"{T.name} harness_flags")
    opt = tval(T, cfg, "harness_optimize", "-O2")
    base = (f"{ctx['clang']} -std={tval(T, cfg, 'std', 'c++14')} {opt} -ffp-contract=off -fplugin={ctx['plugin']} "
            f"-mllvm -enzyme-loose-types=1 -mllvm -enzyme-globals-default-inactive=1 "
            f"-Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas -Wno-deprecated-declarations "
            f"{defs} -I{W / 'gen'} -I{ROOT / 'lib'} -I{T.dir} {T.dir / tval(T, cfg, 'harness', 'harness.cc')} "
            f"{flags} -o {W / 'harness_bin'}")
    rc = sh(S, base, W / "compile.log", timeout=1800, harness=True)
    fwd = False
    if rc:
        print(f"compile {cfg}: full build failed (exit {rc}, {W / 'compile.log'}); retrying forward-only: reverse rows will be REVERSE SKIPPED")
        rc = sh(S, base + " -DFWD_ONLY", W / "compile_fwd.log", timeout=1800, harness=True)
        fwd = True
        if rc:
            print(tail(W / "compile_fwd.log", 15))
            raise Fail(f"harness compile of {cfg} failed", 4)
    args = subst(tval(T, cfg, "harness_args", "{cfg}"), ctx, f"{T.name} harness_args")
    bars = T.dir / tval(T, cfg, "bars", "bars.txt")
    rc = sh(S, f"COFACTOR_BARS={bars} {W / 'harness_bin'} {args}", W / "results.txt",
            timeout=int(tval(T, cfg, "harness_timeout", "3600")), harness=True)
    res = Path(W / "results.txt").read_text().split("\n")
    summ = [l for l in res if l.startswith("cofactor ") or "| FAIL" in l or "BAR WITHOUT ROW" in l]
    print("\n".join(summ) if summ else tail(W / "results.txt"))
    print(f"harness {cfg}: {WORD.get(rc, f'exit {rc}')}{' (forward-only build)' if fwd else ''}, results {W / 'results.txt'}")
    return rc if rc in (0, 1, 2) else 4


def do_run(S, T, cfgs, src, accept=False):
    require_batch(S, "run")
    res = {}
    for cfg in cfgs:
        print(f"===== {T.name} / {cfg}")
        try:
            s = resolve_src(T, cfg, src, need=True)
            do_recipe(S, T, cfg, s, "build", "build", "build")
            do_recipe(S, T, cfg, s, "dump", "dumprun", "dump")
            res[cfg] = do_harness(S, T, cfg, s, accept)
        except Fail as e:
            print(f"cofactor: {e}")
            res[cfg] = e.code
    print("===== summary")
    exp = {c: tval(T, c, "expect", "PASS") for c in cfgs}
    for cfg, rc in res.items():
        note = "" if WORD.get(rc) == exp[cfg] else f"   <-- expected {exp[cfg]}"
        print(f"  {T.name} {cfg:24s} {WORD.get(rc, rc)}{note}")
    bad = [c for c, rc in res.items() if WORD.get(rc) != exp[c]]
    print(f"cofactor run {T.name}: {len(res) - len(bad)}/{len(res)} configs as expected")
    return 0 if not bad else max(res[c] for c in bad) or 1


def job_banner(commit, dirty, frozen):
    """the job script's first line: which cofactor runs, shell-quoted (commit messages may contain quotes)"""
    return "echo " + shlex.quote(f"cofactor {commit}{' + uncommitted changes' if dirty else ''} (frozen copy {frozen})")


def do_submit(S, T, cfgs, src, accept=False):
    base = Path(os.path.expandvars(S.get("paths", "work"))) / T.name
    base.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    # the job runs frozen copies: later edits of cofactor or of the project's .cofactor/ cannot reach a queued job
    frozen = base / f"cofactor_{stamp}"
    shutil.copytree(ROOT, frozen, ignore=shutil.ignore_patterns(".git", "__pycache__", "site.ini"))
    if T.project is not None:
        fp = base / f"frozen_{stamp}" / T.project.name  # keeps the project's name: same work and manifest directories as unfrozen runs
        shutil.copytree(T.project / ".cofactor", fp / ".cofactor", ignore=shutil.ignore_patterns("__pycache__"))
        tdir = fp / T.dir.relative_to(T.project)
    else:
        tdir = base / f"frozen_{stamp}" / T.dir.name
        shutil.copytree(T.dir, tdir, ignore=shutil.ignore_patterns("__pycache__"))
    job = base / f"job_{stamp}.sh"
    src = resolve_src(T, cfgs[0], src, need=True)
    cmd = [sys.executable, str(frozen / "bin" / "cofactor"), "run", str(tdir), "-s", str(src)]
    for c in cfgs:
        cmd += ["-c", c]
    if accept:
        cmd.append("--accept")
    hdr = "\n".join(lines(S.get("batch", "header", fallback="")))
    commit = subprocess.run(["git", "-C", str(ROOT), "log", "--oneline", "-1"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--short"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    job.write_text(f"#!/bin/bash\n#SBATCH -J cof_{T.name.replace('/', '_')}\n#SBATCH -o {base}/job_{stamp}_%j.out\n{hdr}\n"
                   f"{job_banner(commit, bool(dirty), frozen)}\n"
                   f"export COFACTOR_SITE={S.path}\n{' '.join(shlex.quote(x) for x in cmd)}\n")
    r = subprocess.run(S.get("batch", "submit", fallback="sbatch").split() + [str(job)], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    print(r.stdout.strip())
    print(f"job script {job}; output {base}/job_{stamp}_<jobid>.out")
    return r.returncode


def do_check_env(S):
    print(f"site file: {S.path}")
    plugin = Path(os.path.expandvars(S.get("toolchain", "enzyme_plugin", fallback="")))
    print(f"enzyme plugin: {plugin} {'found' if plugin.is_file() else 'MISSING'}")
    work = Path(os.path.expandvars(S.get("paths", "work", fallback="")))
    print(f"work directory: {work} {'exists' if work.is_dir() else '(created on first use)'}")
    rb = S.get("batch", "require_batch", fallback="no")
    print(f"require_batch: {rb}; inside a batch job now: {bool(os.environ.get(S.get('batch', 'job_env', fallback='SLURM_JOB_ID')))}")
    P = find_project(START)
    print(f"project: {P if P else '(none: not inside a directory tree with .cofactor/)'}")
    print(f"checks: {', '.join(all_targets())}")
    print("smoke test: 'cofactor submit examples/c_cooling/cooling' (or 'cofactor run ...' where computing is allowed)")
    return 0 if plugin.is_file() else 4


def do_list():
    P = find_project(START)
    print(f"project: {P if P else '(none)'}")
    for name, d in all_targets().items():
        T = load_target(str(d))
        print(f"{name}: {T.get('target', 'description', fallback='')}")
        print(f"    configs: {', '.join(configs(T))}")
    return 0


def do_init(where):
    d = Path(where).resolve() / ".cofactor"
    if d.exists():
        raise Fail(f"{d} exists")
    (d / "checks").mkdir(parents=True)
    shutil.copy(ROOT / "templates" / "project.ini", d / "project.ini")
    shutil.copytree(ROOT / "templates" / "check", d / "checks" / "example")
    print(f"created {d}: project.ini and checks/example (docs/WRITING_A_CHECK.md); commit it with the project")
    return 0


def do_new(name):
    P = find_project(START)
    if P is None:
        raise Fail("not inside a project with .cofactor/: run 'cofactor init' at the project root first")
    d = P / ".cofactor" / "checks" / name
    if d.exists():
        raise Fail(f"{d} exists")
    shutil.copytree(ROOT / "templates" / "check", d)
    print(f"created {d}: edit target.ini, spec.py, harness.cc, bars.txt (docs/WRITING_A_CHECK.md)")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cofactor", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    p = sub.add_parser("init")
    p.add_argument("-C", dest="where", default=".", help="project root (default: current directory)")
    p = sub.add_parser("new")
    p.add_argument("name")
    sub.add_parser("check-env")
    sub.add_parser("list")
    for c in ("extract", "build", "dump", "harness", "run", "submit"):
        p = sub.add_parser(c)
        p.add_argument("target")
        p.add_argument("-s", "--src", help="source tree of the code under test (default: the target's default_src)")
        p.add_argument("-c", "--config", action="append", help="config (repeatable; run/submit default: all)")
        p.add_argument("--accept", action="store_true", help="accept a changed extraction (rewrite the manifest)")
    p = sub.add_parser("compare")
    p.add_argument("old")
    p.add_argument("new")
    a = ap.parse_args(argv)
    try:
        if a.cmd is None:
            ap.print_help()
            return 0
        if a.cmd == "list":
            return do_list()
        if a.cmd == "init":
            return do_init(a.where)
        if a.cmd == "new":
            return do_new(a.name)
        if a.cmd == "compare":
            return subprocess.call([sys.executable, str(ROOT / "tools" / "compare_results.py"), a.old, a.new])
        S = load_site()
        if a.cmd == "check-env":
            return do_check_env(S)
        T = load_target(a.target)
        cfgs = a.config or configs(T)
        bad = [c for c in cfgs if c not in configs(T)]
        if bad:
            raise Fail(f"unknown config(s) {bad}; target {T.name} has {configs(T)}")
        if a.cmd in ("extract", "build", "dump", "harness") and len(cfgs) != 1:
            raise Fail(f"'{a.cmd}' takes exactly one -c CONFIG")
        if a.cmd == "submit":
            resolve_src(T, cfgs[0], a.src, need=True)
            return do_submit(S, T, cfgs, a.src, a.accept)
        if a.cmd == "run":
            return do_run(S, T, cfgs, a.src, a.accept)
        src = resolve_src(T, cfgs[0], a.src, need=True)
        if a.cmd == "extract":
            do_extract(S, T, cfgs[0], src, a.accept)
        elif a.cmd == "build":
            do_recipe(S, T, cfgs[0], src, "build", "build", "build")
        elif a.cmd == "dump":
            do_recipe(S, T, cfgs[0], src, "dump", "dumprun", "dump")
        elif a.cmd == "harness":
            return do_harness(S, T, cfgs[0], src, a.accept)
        return 0
    except Fail as e:
        print(f"cofactor: {e}")
        return e.code
