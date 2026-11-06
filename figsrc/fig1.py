"""Fig. 1 (single column, page 1): three ways to hold many rhythms.

Left: standalone clocks, one-way rings that share one state, a ring per rhythm.
Middle: one ring, whose rhythms are overtones of its single cycle.
Right: nested cycles, one chain of ever slower states with shortcuts back to the start, a shortcut per rhythm.
Every cycle passes through the gray state; darker states are slower. Each sketch has ten states."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from style import *

HERE = os.path.dirname(os.path.abspath(__file__))

X0, X1, Y0, Y1 = -0.90, 4.78, -1.40, 1.00
W = COLW
H = round(W * (Y1 - Y0) / (X1 - X0), 3)
fig = figure(W, H)
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_axis_off()
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect("equal")

RN, RG = 0.075, 0.095          # radius of a state, of the gray state
SHARED = ("#CFCFCF", GRAY)
YT, YS, LS = -0.90, -1.10, 0.185   # names, first line of the descriptions, line step


def node(x, y, r, fc, ec, lw=0.7, z=3):
    ax.add_patch(Circle((x, y), r, fc=fc, ec=ec, lw=lw, zorder=z))


def arrow(p, q, col, lw=1.1, ms=7.0, shrinkA=5.0, shrinkB=5.4, z=2, conn="arc3"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=col, lw=lw,
                                 shrinkA=shrinkA, shrinkB=shrinkB, zorder=z, connectionstyle=conn))


def arc_arrow(c, R, a0, a1, r0, r1, col, ms=6.5, lw=1.1, gap=0.04):
    """arrow along the circle (center c, radius R) from the state at angle a0 to the state at angle a1, with equal
    clearance from both state outlines; the quadratic Bezier of arc3 follows the circle (sagitta matched)"""
    sgn = np.sign(a1 - a0)
    b0, b1 = a0 + sgn * (r0 + gap) / R, a1 - sgn * (r1 + gap) / R
    p = (c[0] + R * np.cos(b0), c[1] + R * np.sin(b0))
    q = (c[0] + R * np.cos(b1), c[1] + R * np.sin(b1))
    half = abs(b1 - b0) / 2
    rad = sgn * (1 - np.cos(half)) / np.sin(half)              # 2 * sagitta / chord
    arrow(p, q, col, lw=lw, ms=ms, shrinkA=0, shrinkB=0, conn=f"arc3,rad={rad:.5f}")


def names(x, name, what, col):
    ax.text(x, YT, name, color=col, fontsize=8, ha="center", va="center")
    for i, line in enumerate(what):
        ax.text(x, YS - i * LS, line, color=GRAY, fontsize=7, ha="center", va="center")


# ---------------------------------------------------------------- standalone clocks: three rings share one state
rr = 0.42
tones = [("#A9CCF5", RING), (RING_LIGHT, RING), (RING, RING_DARK)]     # fast -> slow rings
for (fc, ec), ang in zip(tones, (90, 210, 330)):
    th = np.radians(ang)
    cx, cy = rr * np.cos(th), rr * np.sin(th)
    phis = th + np.pi + np.arange(4) * 2 * np.pi / 4
    pts = [(cx + rr * np.cos(p), cy + rr * np.sin(p)) for p in phis]
    for k in range(4):
        arrow(pts[k], pts[(k + 1) % 4], RING, conn="arc3,rad=-0.25")
    for k in range(1, 4):
        node(*pts[k], RN, fc, ec)
node(0, 0, RG, *SHARED)
names(0.0, "standalone clocks", ("one rhythm", "per ring"), RING)

# ---------------------------------------------------------------- one ring: equal rates, rhythms are overtones
R1, cx1, cy1, n1 = 0.56, 1.705, 0.13, 10
phis = np.pi / 2 + np.arange(n1) * 2 * np.pi / n1                      # gray state on top, one-way around
pts = [(cx1 + R1 * np.cos(p), cy1 + R1 * np.sin(p)) for p in phis]
for k in range(n1):                                                     # each arrow midway between its states
    arc_arrow((cx1, cy1), R1, phis[k], phis[k] + 2 * np.pi / n1, RG if k == 0 else RN, RG if k == n1 - 1 else RN,
              ONE, ms=6.5)
for k in range(1, n1):
    node(*pts[k], RN, ONE_LIGHT, ONE)
node(*pts[0], RG, *SHARED)
names(cx1, "one ring", ("overtones of", "one cycle"), ONE)

# ---------------------------------------------------------------- nested cycles: ever slower states, shortcuts home
x0, dx, s, nb = 2.705, 0.215, 3, 3
xs = [x0 + dx * j for j in range(s * nb + 1)]
yy = -0.15
oranges = plt.matplotlib.colors.LinearSegmentedColormap.from_list("o", [NEST_LIGHT, NEST, NEST_DARK])  # fast -> slow
for j in range(1, len(xs)):
    arrow((xs[j - 1], yy), (xs[j], yy), NEST, ms=6.0, shrinkA=4.4, shrinkB=4.7)
for blk in range(1, nb + 1):
    j = blk * s
    arrow((xs[j], yy + 0.07), (xs[0] + 0.03, yy + 0.08), NEST, shrinkA=1.5, shrinkB=4.8,
          conn=f"arc3,rad={0.42 + 0.03 * blk}")
for j in range(1, len(xs)):
    node(xs[j], yy, RN, oranges((j - 1) / (len(xs) - 2)), NEST_DARK)
node(xs[0], yy, RG, *SHARED)
names((xs[0] + xs[-1]) / 2, "nested cycles", ("one rhythm", "per shortcut"), NEST)

save(fig, os.path.join(HERE, "out", "fig1.pdf"))
print("fig1 size in", round(W, 3), H)
