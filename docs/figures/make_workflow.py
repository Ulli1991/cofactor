"""Draws docs/figures/workflow.png: the cofactor workflow and how a bug shows up (GIZMO C6 softened gravity, issue #178).
usage: python3 make_workflow.py <out.png>"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

W, H = 980, 560
fig = plt.figure(figsize=(W / 100, H / 100), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")

C = {"box": ("#f3f6fa", "#4a6a8a"), "hl": ("#fff4e5", "#c07a1a"), "bad": ("#fdecec", "#b03a3a"),
     "ok": ("#eaf6ec", "#3a8a4a"), "opt": ("#ffffff", "#4a6a8a")}
INK, SUB = "#1b2a38", "#3c4c5c"


def box(x, y, w, h, kind="box"):
    fc, ec = C[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=8", fc=fc, ec=ec, lw=1.3,
                                ls="--" if kind == "opt" else "-"))


def text(x, y, s, size=9.0, bold=False, mono=False, color=None):
    ax.text(x, y, s, fontsize=size, fontweight="bold" if bold else "normal", va="baseline",
            family="monospace" if mono else "sans-serif", color=color or (INK if (bold or mono) else SUB))


def arrow(x0, y0, x1, y1, dashed=False):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11, lw=1.3, color="#4a6a8a",
                                 ls="--" if dashed else "-"))


text(20, 30, "How cofactor checks a physics implementation", 12, True)

box(208, 48, 364, 58, "opt")
text(222, 70, "optional: states from a real run", 9.6, True)
text(222, 90, "test-only hook in the code  →  build (no FMA)  →  dump run", 8.2)
arrow(440, 106, 440, 138, dashed=True)

steps = [
    (20, 168, "box", "1  Your code", ["the real source tree", "e.g. potential φ(u) and", "its hand-coded force", "", "checks live in .cofactor/"]),
    (216, 168, "box", "2  Extract", ["copy functions verbatim", "by unique markers", "hash every block: stops", "when the source moved"]),
    (412, 168, "box", "3  Harness", ["compile with clang +", "Enzyme (no FMA)", "row 0: the copy", "reproduces the code", "bit for bit"]),
    (608, 168, "hl", "4  Compare", ["the code's claims vs", "Enzyme on that code:", "derivatives, Jacobians,", "adjoints, identities,", "energy bookings"]),
    (804, 156, "box", "5  Report", ["one row per check:", "worst value + where", "PASS / FAIL against", "bars set beforehand"]),
]
for x, w, kind, head, lines in steps:
    box(x, 140, w, 126, kind)
    text(x + 14, 164, head, 9.8, True)
    for i, l in enumerate(lines):
        text(x + 14, 186 + 16.5 * i, l, 8.2)
for x in (188, 384, 580, 776):
    arrow(x, 203, x + 26, 203)
box(804, 284, 156, 40, "opt")
text(818, 309, "bars.txt (fixed first)", 8.2)
arrow(882, 284, 882, 268, dashed=True)

text(20, 366, "How a bug shows up: GIZMO Wendland C6 softened gravity (issue #178)", 12, True)
panels = [
    (20, 220, "box", "The code claims", [("mode 1 = (1/u) dφ/du", True)],
     ["force (mode 1) and", "potential φ (mode −1) are two", "separate hand-written", "polynomials in kernel.h"]),
    (262, 220, "bad", "cofactor row (g1)", [("max_rel = 0.75", True), ("bar 1e-12  →  FAIL", True)],
     ["Enzyme: (1/u) dφ/du is exactly", "4 × the coded force; and", "φ(1⁻) = −4, not −1 (row g5)"]),
    (504, 220, "hl", "Rule out the harness", [],
     ["force passes Poisson vs the", "SPH kernel, only φ is off", "control arm with 1/64", "instead of 1/16 passes", "reproducer: softening 4×", "the requested value"]),
    (746, 214, "ok", "Fix", [("- (1./16.) * (...)", True), ("+ (1./64.) * (...)", True)],
     ["φ(1) = −1, force =", "(1/u) dφ/du to round-off,", "softening as requested"]),
]
for x, w, kind, head, mono, lines in panels:
    box(x, 384, w, 150, kind)
    text(x + 14, 408, head, 9.8, True)
    y = 432
    for s, _ in mono:
        text(x + 14, y, s, 8.2, mono=True)
        y += 18
    y += 6 if mono else 0
    for l in lines:
        text(x + 14, y, l, 8.2)
        y += 16.5
for x in (240, 482, 724):
    arrow(x, 459, x + 20, 459)

fig.savefig(sys.argv[1], dpi=200, facecolor="white")
print("wrote", sys.argv[1])
