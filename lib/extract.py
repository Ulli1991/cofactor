"""cofactor: copy real GIZMO code into a harness, verbatim, by text markers, and fail loudly when the source moved.

Rules (every one is an error, exit 3):
- a start marker must occur exactly once in its search range, a stop marker at least once after it;
- a function ends at its matching closing brace (strings, chars and comments skipped), not at the first '}' in column 0;
- preprocessor repairs (#endif dropped / #if closed) must equal what the spec declares;
- every extracted block is hashed into EXTRACTED.txt; if the target keeps a manifest (manifest_<config>.txt) and a hash differs,
  the run stops with "target changed" unless --accept is given (then the manifest is rewritten and the change logged).
"""
import hashlib
import re
import sys
from pathlib import Path


class ExtractError(Exception):
    pass


def brace_depth_change(lines):
    """net { minus } over the lines, skipping // and /* */ comments, string and char literals"""
    depth, in_c = 0, False
    for l in lines:
        i, n = 0, len(l)
        while i < n:
            c = l[i]
            if in_c:
                if l.startswith("*/", i):
                    in_c, i = False, i + 2
                    continue
                i += 1
                continue
            if l.startswith("//", i):
                break
            if l.startswith("/*", i):
                in_c, i = True, i + 2
                continue
            if c in "\"'":
                q, i = c, i + 1
                while i < n and l[i] != q:
                    i += 2 if l[i] == "\\" else 1
                i += 1
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
            i += 1
    return depth


class Extractor:
    def __init__(self, build_dir, out_dir, target_dir=None):
        self.B, self.O = Path(build_dir), Path(out_dir)
        self.T = Path(target_dir) if target_dir else None
        self.O.mkdir(parents=True, exist_ok=True)
        self.rows = []  # (name, where, sha256)

    def _lines(self, rel):
        p = self.B / rel
        if not p.is_file():
            raise ExtractError(f"{rel}: file missing in {self.B}")
        return p.read_text().split("\n")

    @staticmethod
    def _find_once(L, text, lo, hi, what):
        hits = [i for i in range(lo, hi) if text in L[i]]
        if len(hits) != 1:
            raise ExtractError(f"{what} marker {text.strip()[:70]!r}: {len(hits)} matches in lines {lo + 1}-{hi} (need exactly 1)")
        return hits[0]

    def _balance(self, L, expect_dropped, expect_closed, where):
        out, depth, dropped = [], 0, 0
        for l in L:
            s = l.strip()
            if s.startswith("#if"):
                depth += 1
            elif s.startswith("#endif"):
                if depth == 0:
                    dropped += 1
                    continue
                depth -= 1
            out.append(l)
        if (dropped, depth) != (expect_dropped, expect_closed):
            raise ExtractError(f"{where}: preprocessor repair dropped {dropped} #endif / closed {depth} #if; "
                               f"spec expects {expect_dropped} / {expect_closed}")
        return out + ["#endif"] * depth

    def _emit(self, name, where, text):
        (self.O / name).write_text(text)
        self.rows.append((name, where, hashlib.sha256(text.encode()).hexdigest()[:16]))

    def block(self, name, rel, start, stop, include_stop=False, start_after=None, expect_dropped=0, expect_closed=0,
              cut_last=None, prefix="", suffix=""):
        """lines from the unique start marker up to (excluding, or including) the first stop marker after it;
        cut_last: then drop everything from the last occurrence of this text (must occur); prefix/suffix: literal wrap"""
        L = self._lines(rel)
        lo = 0
        if start_after is not None:
            lo = self._find_once(L, start_after, 0, len(L), f"{rel}: start_after") + 1
        a = self._find_once(L, start, lo, len(L), f"{rel}: start")
        b = next((i for i in range(a + 1, len(L)) if stop in L[i]), None)
        if b is None:
            raise ExtractError(f"{rel}: stop marker {stop.strip()[:70]!r} not found after line {a + 1}")
        if include_stop:
            b += 1
        where = f"{rel}:{a + 1}-{b}"
        text = "\n".join(self._balance(L[a:b], expect_dropped, expect_closed, where)) + "\n"
        if cut_last is not None:
            if cut_last not in text:
                raise ExtractError(f"{where}: cut_last marker {cut_last!r} not in block")
            text = text.rsplit(cut_last, 1)[0]
        self._emit(name, where, prefix + text + suffix)

    def function(self, rel, sig_regex, newsig, nth=0):
        """a top-level function: signature line replaced by newsig, body verbatim to its matching brace; returns the text"""
        L = self._lines(rel)
        hits = [i for i, l in enumerate(L) if re.match(sig_regex, l)]
        if nth >= len(hits):
            raise ExtractError(f"{rel}: signature {sig_regex!r}: {len(hits)} matches, need index {nth}")
        a = hits[nth]
        open_on_sig = brace_depth_change([L[a]]) > 0
        if not open_on_sig and not L[a + 1].lstrip().startswith("{"):
            raise ExtractError(f"{rel}:{a + 1}: signature not followed by its body")
        depth, b = 0, None
        for i in range(a, len(L)):
            depth += brace_depth_change([L[i]])
            if i > a and depth == 0:
                b = i
                break
        if b is None:
            raise ExtractError(f"{rel}:{a + 1}: no matching closing brace for {sig_regex!r}")
        body = L[a + 1:b + 1]
        if "{arg}" in newsig:  # keep the original name of the single int argument
            m = re.search(r"\(int (\w+)\)", L[a])
            if m is None:
                raise ExtractError(f"{rel}:{a + 1}: newsig uses {{arg}} but the signature has no single int argument")
            newsig = newsig.replace("{arg}", m.group(1))
        return f"{rel}:{a + 1}-{b + 1}", newsig + (" {" if open_on_sig else "") + "\n" + "\n".join(body) + "\n"

    def functions(self, name, items):
        """several function() results into one file; items: (rel, sig_regex, newsig[, nth]) or literal str lines"""
        parts, wh = [], []
        for it in items:
            if isinstance(it, str):
                parts.append(it)
                continue
            where, text = self.function(*it)
            parts.append(text)
            wh.append(where)
        self._emit(name, ", ".join(wh) or "(literal)", "\n".join(parts))

    def literal(self, name, text):
        self._emit(name, "(literal)", text)

    def finish(self, accept=False):
        man = "\n".join(f"{n}  {h}  {w}" for n, w, h in self.rows) + "\n"
        (self.O / "EXTRACTED.txt").write_text(man)
        if self.T is None:
            return 0
        self.T.mkdir(parents=True, exist_ok=True)
        mp = self.T / f"manifest_{self.O.name}.txt"  # one manifest per build config (out_dir name)
        if not mp.is_file():
            mp.write_text(man)
            print(f"cofactor: manifest created ({len(self.rows)} blocks)")
            return 0
        old = {l.split()[0]: l.split()[1] for l in mp.read_text().split("\n") if l.strip()}
        new = {n: h for n, w, h in self.rows}
        changed = sorted(k for k in set(old) | set(new) if old.get(k) != new.get(k))
        if not changed:
            print(f"cofactor: extraction unchanged ({len(self.rows)} blocks)")
            return 0
        print("cofactor: target changed: " + ", ".join(changed))
        if accept:
            mp.write_text(man)
            print("cofactor: manifest rewritten (--accept)")
            return 0
        return 3


def run(spec, argv):
    """spec(X, build_dir) fills the extractor; argv: build_dir out_dir target_dir [--accept]"""
    if len(argv) < 4:
        print("usage: spec.py <build_dir> <out_dir> <target_dir> [--accept]")
        sys.exit(2)
    X = Extractor(argv[1], argv[2], argv[3])
    try:
        spec(X, Path(argv[1]))
    except (ExtractError, StopIteration) as e:
        print(f"cofactor: EXTRACTION ERROR: {e}")
        sys.exit(3)
    for n, w, h in X.rows:
        print(f"  {n:24s} {h}  {w}")
    sys.exit(X.finish(accept="--accept" in argv))
