"""Fig. 4 (single column, two panels): rhythms cost states and drive, and measured rhythms reveal both.

(a) States and drive per rhythm for Q > 10, relative to standalone clocks: both shared-budget bounds settle at one
    half; the single-rhythm drive bound falls off like 1/N, the single-rhythm state bound levels off at
    2/(n0 - 1), about 0.065; standalone clocks sit at one.
(b) Hearing a hidden clock: cross-correlation of a three-valued readout, its fit, and its three rhythms."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from style import *

HERE = os.path.dirname(os.path.abspath(__file__))
ab = np.load(os.path.join(HERE, "data", "fig4a_costs.npz"))
tr = np.load(os.path.join(HERE, "data", "fig4b_trace.npz"))

W = COLW
L, R, GAPX = 0.34, 0.03, 0.48
PW = (W - L - R - GAPX) / 2
B, PH, TOP = 0.29, 1.55, 0.17
Hh = B + PH + TOP
fig = figure(W, Hh)
axa = axes_in(fig, L, B, PW, PH)
hb2 = 0.44
axb2 = axes_in(fig, L + PW + GAPX, B, PW, hb2)                        # the three rhythms
axb = axes_in(fig, L + PW + GAPX, B + hb2 + 0.03, PW, PH - hb2 - 0.03)  # data and fit
XR = L + PW + GAPX - 0.38

# ------------------------------------------------------------------ (a) cost per rhythm, relative to clocks
Q0 = 10.0
N = ab["N"].astype(float)
n0 = int(ab["n0_inf"])                                  # ring size of one standalone clock (32 for Q > 10)
clock_states = N * (n0 - 1) + 1                          # N standalone clocks sharing one state
states = ab["multi_inf"] * N / clock_states             # capacity, Eq. (2) in its sharp form
states_single = ab["single_inf"] * N / clock_states     # d >= max(n0, 2N + 1)
drive = ab["aff_multi"] / (2 * np.pi * Q0)              # drive law, discrete form
drive_single = ab["aff_single"] / (2 * np.pi * Q0)      # 1/N: one ring needs 2 pi Q
single_lo = np.minimum(states_single, drive_single)
single_hi = np.maximum(states_single, drive_single)
shared = np.minimum(np.maximum(states, states_single), np.maximum(drive, drive_single))
ax = axa
ax.fill_between(N, 0, single_lo, color=LGRAY, lw=0, zorder=0)
ax.fill_between(N, single_hi, shared, where=shared > single_hi, facecolor="none", edgecolor=MGRAY,
                hatch="//////", lw=0, zorder=0.5, interpolate=True)
ax.axhline(1.0, color=RING, lw=1.3, ls=(0, (4.5, 2)), zorder=2)
ax.plot(N, states_single, ls=(0, (1, 1.3)), color=GRAY, lw=1.1, zorder=2)
ax.plot(N, drive_single, ls=(0, (1, 1.3)), color=GRAY, lw=1.1, zorder=2)
ax.plot(N, drive, color=TURN, lw=1.3, zorder=3)
ax.plot(N, states, color=INK, lw=0.7, zorder=4)
ax.plot(N, states, ls="none", marker="o", ms=2.4, mfc=INK, mec=INK, mew=0.5, zorder=5)
ax.set_xlim(0, 31); ax.set_ylim(0, 1.17)
ax.set_xticks([0, 10, 20, 30]); ax.set_yticks([0, 0.5, 1.0]); ax.set_yticklabels(["0", "0.5", "1"])
minor(ax)
ax.text(16.0, 0.955, "standalone clocks", color=RING, fontsize=7, ha="center", va="top")
ax.text(9.0, 0.635, "drive", color=TURN, fontsize=7, ha="center", va="bottom")
ax.text(16.0, 0.445, "states", color=INK, fontsize=7, ha="center", va="top",
        bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
ax.text(11.5, 0.195, "single rhythm", color=GRAY, fontsize=7, ha="left", va="center",
        bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
ax.set_xlabel(r"Rhythms $N$")
ax.set_ylabel(r"Relative cost", labelpad=1.5)
header(fig, ax, "(a)", r"Cost per rhythm, $Q>10$", dy=0.055)

# ------------------------------------------------------------------ (b) hearing a hidden clock
tau = tr["tau"]; C = tr["C"]; lo = tr["lo"]; hi = tr["hi"]
tt = np.linspace(0, 51.24, 1500)
comps = [cc[0] * np.exp(-cc[1] * tt) * np.cos(cc[2] * tt + cc[3]) for cc in tr["comps"]]
fit = sum(comps) + sum(bb * np.exp(-cc * tt) for bb, cc in tr["exps"])
ax = axb
ax.axhline(0, color=MGRAY, lw=0.5, zorder=0)
ax.vlines(tau, lo, hi, color=DATA, lw=0.7, zorder=2)
ax.plot(tau, C, ls="none", marker="o", ms=2.5, mfc="white", mec=DATA, mew=0.7, zorder=3)
ax.plot(tt, fit, color=INK, lw=1.0, zorder=4)
ax.set_xlim(-2, 56); ax.set_ylim(-0.135, 0.068)
ax.set_yticks([-0.1, 0.0]); ax.set_yticklabels([r"$-0.1$", r"$0$"])
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
ax.set_xticks([0, 20, 40]); ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
ax.tick_params(labelbottom=False)
ax.text(54, 0.050, "data and fit", fontsize=7, ha="right", va="top")
ax.set_ylabel(r"$C_{10}(\tau)$", labelpad=1.0)
ax.yaxis.set_label_coords(-0.33 / PW, 0.42)              # midway between panel (a) and the tick labels
header(fig, ax, "(b)", "Hearing a hidden clock", dy=0.055)
ax = axb2
ax.axhline(0, color=MGRAY, lw=0.5, zorder=0)
for cc, col in zip(comps, SEQ_RING):
    ax.plot(tt, cc, color=col, lw=1.1)
ax.set_xlim(-2, 56); ax.set_ylim(-0.05, 0.05)
ax.set_yticks([]); ax.set_xticks([0, 20, 40]); ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
ax.tick_params(axis="y", which="both", left=False, right=False)
ax.tick_params(axis="x", which="both", top=False)
ax.text(54, 0.041, "three rhythms", fontsize=7, ha="right", va="top", color=RING)
ax.set_xlabel(r"Lag $\tau$")

save(fig, os.path.join(HERE, "out", "fig4.pdf"))
print("fig4 size in", round(W, 3), round(Hh, 3))
