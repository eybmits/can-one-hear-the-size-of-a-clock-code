"""Fig. 3 (single column, two panels): sharing halves the state cost at most, and only by blurring.

(a) Weyl law for clocks: the most rhythms sharper than Q per state, times Q, with standalone clocks (blue), one
    large ring (pink: arccot(Q)/pi rhythms per state, exactly half the capacity for Q > 1, touching the blue
    sawtooth at its lower corners), and nested cycles (orange).
(b) The price of a dip: state cost of nested cycles relative to standalone clocks against the spacing of
    neighboring rhythms in mean half-widths (s = 4 and 32)."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter, FixedLocator
from style import *

HERE = os.path.dirname(os.path.abspath(__file__))
b = np.load(os.path.join(HERE, "data", "fig3a_density.npz"))
c = np.load(os.path.join(HERE, "data", "fig3b_spacing.npz"))

W = COLW
L, R, GAPX = 0.34, 0.03, 0.35
PW = (W - L - R - GAPX) / 2
B, PH, TOP = 0.29, 1.55, 0.17
Hh = B + PH + TOP
fig = figure(W, Hh)
axa = axes_in(fig, L, B, PW, PH)
axb = axes_in(fig, L + PW + GAPX, B, PW, PH)
XR = L + PW + GAPX - 0.35

# ------------------------------------------------------------------ (a) rhythm density
ax = axa
Q0 = b["Q0"]
ax.fill_between(Q0, b["bound"], 0.8, color=LGRAY, lw=0, zorder=0)
ax.axhline(2 / np.pi, color=GRAY, lw=0.5, ls=(0, (1, 1.2)), zorder=1)
ax.axhline(1 / np.pi, color=GRAY, lw=0.5, ls=(0, (1, 1.2)), zorder=1)
ax.plot(Q0, Q0 * np.arctan(1 / Q0) / np.pi, color=ONE, lw=1.2, zorder=2.5)   # one large ring
ax.plot(Q0, b["rings"], color=RING, lw=1.2, zorder=3)
ax.plot(Q0, b["chain"], color=NEST, lw=1.2, zorder=4)
ax.plot(Q0, b["bound"], color=INK, lw=1.2, zorder=5)
ax.set_xscale("log")
ax.set_xlim(0.3, 40); ax.set_ylim(0, 0.8)
ax.set_xticks([1, 10]); ax.set_xticklabels(["1", "10"])
ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8]); ax.set_yticklabels(["0", "", "0.4", "", "0.8"])
minor(ax, x=False)
ax.text(33, 2 / np.pi + 0.024, r"$2/\pi$", color=GRAY, fontsize=7, ha="right", va="bottom")
ax.text(33, 1 / np.pi - 0.022, r"$1/\pi$", color=GRAY, fontsize=7, ha="right", va="top")
ax.text(1.25, 0.20, "standalone clocks", color=RING, fontsize=7, ha="left", va="top")
ax.text(0.36, 0.105, "one ring", color=ONE, fontsize=7, ha="left", va="top")
ax.text(2.6, 0.445, "nested cycles", color=NEST, fontsize=7, ha="left", va="top")
ax.text(0.36, 0.662, "capacity", color=INK, fontsize=7, ha="left", va="bottom")
ax.set_xlabel(r"Sharpness threshold $Q$")
ax.set_ylabel(r"$Q\,\nu(Q)$", labelpad=1.5)
header(fig, ax, "(a)", "Weyl law for clocks", dy=0.055)

# ------------------------------------------------------------------ (b) the price of a dip
ax = axb
res = 2 / np.sqrt(3)                                   # spacing at which a dip appears between equal lines
ax.axvspan(res, 2.0, color=DIST_FILL, lw=0, zorder=0)
ax.axhspan(0.3, 0.5, color=LGRAY, lw=0, zorder=0)       # below half a standalone clock: forbidden
ax.axhline(0.5, color=INK, lw=1.1, zorder=2)
ax.axhline(1.0, color=RING, lw=1.2, ls=(0, (4.5, 2)), zorder=2)
YEND = 1.10      # every curve ends 10% above the cost of standalone clocks, after its crossing with them
for s, col in ((4, NEST), (32, NEST_DARK)):
    x, y = c[f"s{s}_x"], c[f"s{s}_y"]
    # Up to YEND both curves rise smoothly. Beyond it, the damping angle for large s starts to oscillate with h
    # (the second Fourier term of G becomes comparable), so the cost stops being a smooth function of the spacing.
    n = int(np.searchsorted(y, YEND))
    xe = x[n - 1] + (YEND - y[n - 1]) * (x[n] - x[n - 1]) / (y[n] - y[n - 1])
    ax.plot(np.r_[x[:n], xe], np.r_[y[:n], YEND], color=col, lw=1.3, zorder=4)
ax.set_xscale("log")
ax.set_xlim(0.035, 1.75); ax.set_ylim(0.35, 1.25)
ax.set_xticks([0.1, 1]); ax.set_xticklabels([r"$0.1$", r"$1$"])
ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([0.5, 1.0]); ax.set_yticklabels(["0.5", "1"])
ax.yaxis.set_minor_locator(FixedLocator([0.4, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2]))
ax.text(0.04, 0.975, "standalone clocks", color=RING, fontsize=7, ha="left", va="top")
ax.text(1.7, 0.53, "distinguishable", color=DIST_TEXT, fontsize=7, ha="right", va="bottom", rotation=90, zorder=6,
        bbox=dict(boxstyle="square,pad=0.05", fc=DIST_FILL, ec="none"))
ax.text(0.118, 0.395, "forbidden", color=FORBID_TEXT, fontsize=7, ha="left", va="baseline")
hand = [Line2D([], [], color=NEST, lw=1.3), Line2D([], [], color=NEST_DARK, lw=1.3)]
ax.legend(hand, [r"$s=4$", r"$s=32$"], loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=2,
          handlelength=1.1, columnspacing=0.7, handletextpad=0.35)
ax.set_xlabel(r"Spacing/half-width $\Delta\omega/\bar a$")
ax.set_ylabel(r"Relative state cost", labelpad=1.5)
header(fig, ax, "(b)", "The price of a dip", dy=0.055)

save(fig, os.path.join(HERE, "out", "fig3.pdf"))
print("fig3 size in", round(W, 3), round(Hh, 3))
