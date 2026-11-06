"""Fig. 2 (full text width, three panels): every tick blurs the clock, so all rhythms of a network share one budget.

(a) A drum: Weyl's law. The l-th Dirichlet mode of a square drum of side pi against l (dots) follows
    omega_l = (2/pi)(1 + sqrt(1 + pi l)) (red, dashed), whose height reveals the area; inset: the mode (3, 2).
(b) A fading rhythm: copies of a one-way ring of 20 states, started together in state 0, drift out of step; the
    return probability P_00(t) oscillates and fades with its envelope e^{-at}. Left window: dials at t = 0, 10, 20,
    40 with the occupation of the states (rim) and the clock hand, the average of e^{2 pi i j/d} over all copies,
    which equals e^{(-a + i omega) t}; circles on the curve mark these times. Right window: the spectral line of the
    sharpest rhythm, at omega with half-width a.
(c) One budget: sharpness Q_k of the k-th sharpest rhythm against k/(d-1), with the capacity, Eq. (2),
    Q_k < 2(d-1)/(pi k), in the red of Weyl's law: standalone clocks (rings of m states, Q = cot(pi/m) at
    k/(d-1) = 1/(m-1)), about a factor of two below (double arrow); a one-way ring of 32 states (on the same line); nested cycles (closest); random networks as a density cloud
    (the regions holding half and 90% of their rhythms).
In-panel text is kept to the essentials; the caption carries the explanations. Style: design study 2026-11-06,
version 1 (refined): envelope band, rounded windows, outlined density cloud, x2 bracket."""
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Wedge
from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter
from scipy.ndimage import gaussian_filter
from style import *

logging.getLogger("matplotlib").setLevel(logging.ERROR)   # equal-aspect dial window: limits adjusted on purpose
HERE = os.path.dirname(os.path.abspath(__file__))
Db = np.load(os.path.join(HERE, "data", "fig2b_ring.npz"))   # one-way ring of 20 states, started in one state
Dc = np.load(os.path.join(HERE, "data", "fig2c_rank.npz"))   # sharpness of the k-th sharpest rhythm

W = TEXTW
B, PH, TOP = 0.30, 1.40, 0.20
b0 = B
Hh = b0 + PH + TOP
fig = figure(W, Hh)

BLUES = LinearSegmentedColormap.from_list("b", ["#EEF3F9", "#D3E5FB", "#5FA2EE", "#1565C0", "#0D3C7A"])
CLOUD = ["#E4E4E4", "#C9C9C9"]                  # random networks: regions holding 90% and 50% of their rhythms
cap = lambda x: 2 / (np.pi * x)                 # noqa: E731  capacity, Eq. (2): Q_k < 2(d-1)/(pi k)
clk = lambda x: 1 / np.tan(np.pi * x / (1 + x))  # noqa: E731  standalone clocks of m = 1 + 1/x states each
d, a1, w1 = int(Db["d"]), float(Db["a1"]), float(Db["w1"])


def arrow(ax, p, q, col, lw=0.6, ms=3.5, style="<->", z=6):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, color=col, lw=lw, shrinkA=0, shrinkB=0,
                                 zorder=z))


def window(axi):
    """A small window inside a panel: a soft rounded frame on a faint background, drawn above the panel."""
    axi.minorticks_off()
    axi.set_zorder(3)
    for sp in axi.spines.values():
        sp.set_visible(False)
    axi.patch.set_alpha(0)
    bb = axi.get_position()
    fig.add_artist(FancyBboxPatch((bb.x0, bb.y0), bb.width, bb.height, boxstyle="round,pad=0,rounding_size=0.006",
                                  transform=fig.transFigure, fc="#FBFCFE", ec="#B9C0C9", lw=0.55, zorder=2))


# ================================================================ (a) a drum: Weyl's law
ax = axes_in(fig, 0.40, b0, 1.35, PH)
mm = np.arange(1, 60)
M, N = np.meshgrid(mm, mm)
om_n = np.sort(np.sqrt(M ** 2 + N ** 2).ravel())[:400]       # Dirichlet modes of a square drum of side pi
n = np.arange(1, len(om_n) + 1)
nn = np.logspace(0, np.log10(400), 200)
ax.plot(nn, 2 / np.pi * (1 + np.sqrt(1 + np.pi * nn)), color=WEYL, lw=1.0, ls=(0, (3, 2)), zorder=4)   # Weyl's law, drawn over the modes
ax.plot(n, om_n, ls="none", marker="o", ms=1.5, mfc=INK, mec="none", alpha=0.9, zorder=3, clip_on=False)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 400); ax.set_ylim(1, 30)
ax.set_xticks([1, 10, 100]); ax.set_xticklabels([r"$1$", r"$10$", r"$100$"])
ax.set_yticks([1, 10]); ax.set_yticklabels([r"$1$", r"$10$"])
for axis in (ax.xaxis, ax.yaxis):
    axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    axis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r"Mode number $\ell$", labelpad=1.5)
ax.set_ylabel(r"Frequency $\omega_\ell$", labelpad=1.0)
header(fig, ax, "(a)", "A drum: Weyl's law")
bb = ax.get_position()                          # inset: one mode of the drum, with its nodal lines
hs_ = 0.40 * bb.height
axm = fig.add_axes([bb.x0 + 0.07 * bb.width, bb.y0 + 0.54 * bb.height, hs_ * Hh / W, hs_])
g = np.linspace(0, np.pi, 160)
Z = np.abs(np.sin(3 * g)[None, :] * np.sin(2 * g)[:, None]) ** 0.8
axm.imshow(Z, cmap=LinearSegmentedColormap.from_list("drum", ["#FFFFFF", "#9A9A9A", "#2E2E2E"]), origin="lower",
           extent=(0, 1, 0, 1), vmin=0, vmax=1, interpolation="bilinear")
axm.set_xticks([]); axm.set_yticks([])
for sp in axm.spines.values():
    sp.set_linewidth(0.5); sp.set_edgecolor("#9AA3AD")
# key inside the panel, in the empty lower right
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=2.6, mfc=INK, mec="none", label="drum modes"),
                   Line2D([], [], color=WEYL, lw=1.0, ls=(0, (3, 2)), label="Weyl's law")],
          loc="lower right", fontsize=7, frameon=False, handlelength=1.6, handletextpad=0.5, borderaxespad=0.4,
          labelspacing=0.3)

# ================================================================ (b) a fading rhythm, with two windows
ax = axes_in(fig, 2.22, b0, 2.40, PH)
ax.set_zorder(1)
t, P = Db["t"], Db["P"]
YTOP = 0.25                                     # head room for the two windows
env = 2 / d * np.exp(-a1 * t)                   # the rhythm oscillates within 1/d +- (2/d) e^{-at}
ax.axhline(1 / d, color=MGRAY, lw=0.6, ls=(0, (3, 2)), zorder=1)
ax.fill_between(t, np.clip(1 / d - env, 0, None), 1 / d + env, color=RING_FILL, alpha=0.55, lw=0, zorder=0.5)
ax.plot(t, 1 / d - env, color=RING, lw=0.7, ls=(0, (3, 2)), alpha=0.8, zorder=2)
ax.plot(t, 1 / d + env, color=RING, lw=1.0, ls=(0, (3, 2)), zorder=2)
ax.plot(t, P, color=INK, lw=1.3, zorder=3)
times, probs = Db["dial_t"], Db["dial_p"]
for tt in times:
    if np.interp(tt, t, P) <= YTOP:            # the start, P = 1, lies above the axes
        ax.plot([tt], [np.interp(tt, t, P)], ls="none", marker="o", ms=3.0, mfc="white", mec=INK, mew=0.6, zorder=4,
                clip_on=False)
ax.text(11.0, 1 / d + 2 / d * np.exp(-a1 * 11.0) + 0.005, r"$e^{-at}$", color=RING, fontsize=7, ha="left", va="bottom")
ax.set_xlim(0, 80); ax.set_ylim(0, YTOP)
ax.set_xticks([0, 40, 80]); ax.xaxis.set_minor_locator(FixedLocator([20, 60]))
ax.set_yticks([0, 0.05]); ax.set_yticklabels(["0", r"$1/d$"])
ax.yaxis.set_minor_locator(FixedLocator([0.10, 0.15, 0.20, 0.25]))
ax.set_xlabel(r"Time $t$")
ax.set_ylabel(r"Return probability", labelpad=1.5)
header(fig, ax, "(b)", "A fading rhythm")
bb = ax.get_position()


def dial(axw, cx, cy, r, p, tt, vmax=0.30, w=0.30):
    norm = PowerNorm(0.55, vmin=0, vmax=vmax)
    for j in range(d):                          # state j sits at 360 j/d degrees, clockwise from the top
        c = 90 - 360 * j / d
        axw.add_patch(Wedge((cx, cy), r, c - 180 / d, c + 180 / d, width=w * r, fc=BLUES(norm(min(p[j], vmax))),
                            ec="white", lw=0.3, zorder=2))
    ri = (1 - w) * r
    axw.add_patch(Circle((cx, cy), r, fc="none", ec="#C9D3DF", lw=0.35, zorder=3))
    axw.add_patch(Circle((cx, cy), ri, fc="none", ec="#DDE3EA", lw=0.4, zorder=1))
    L = np.exp(-a1 * tt) * ri * 0.93            # the hand: <e^{i theta}> = e^{-a t} e^{i omega t}
    ang = np.pi / 2 - w1 * tt
    axw.plot([cx, cx + L * np.cos(ang)], [cy, cy + L * np.sin(ang)], color=INK, lw=0.95, solid_capstyle="round",
             zorder=4)
    axw.add_patch(Circle((cx, cy), 0.055 * r, fc=INK, ec="none", zorder=5))


# left window: copies drift out of step (dials)
axw = fig.add_axes([bb.x0 + 0.075 * bb.width, bb.y0 + 0.56 * bb.height, 0.585 * bb.width, 0.415 * bb.height])
window(axw)
axw.set_xticks([]); axw.set_yticks([])
n = len(times)
axw.set_xlim(0, n); axw.set_ylim(-0.36, 1.02)
axw.set_aspect("equal", adjustable="datalim", anchor="C")
for i, (tt, p) in enumerate(zip(times, probs)):
    dial(axw, i + 0.5, 0.52, 0.40, p, tt)
    axw.text(i + 0.5, -0.02, (r"$t=%d$" % tt) if i == 0 else (r"$%d$" % tt), fontsize=7, ha="center", va="top",
             color=INK)

# right window: the spectral line of the sharpest rhythm
axi = fig.add_axes([bb.x0 + 0.70 * bb.width, bb.y0 + 0.46 * bb.height, 0.285 * bb.width, 0.515 * bb.height])
window(axi)
Om = np.linspace(0, 2 * w1, 801)
line = (a1 / (a1 ** 2 + (Om - w1) ** 2) + a1 / (a1 ** 2 + (Om + w1) ** 2)) / d
axi.fill_between(Om, 0, line, color=RING_LIGHT, alpha=0.45, lw=0)
axi.plot(Om, line, color=RING, lw=1.1)
hmax = float(line.max())
axi.set_xlim(0, 2 * w1); axi.set_ylim(0, 1.42 * hmax)
axi.set_xticks([w1]); axi.set_xticklabels([r"$\omega$"], fontsize=7)
axi.set_yticks([])
axi.tick_params(axis="both", which="both", top=False, right=False, left=False, length=2, pad=1.5, color="#9AA3AD")
arrow(axi, (w1 - a1, hmax / 2), (w1 + a1, hmax / 2), RING_DARK, lw=0.6, ms=3.5)
axi.text(w1 + a1 + 0.02, hmax / 2, r"$2a$", color=RING_DARK, fontsize=7, ha="left", va="center")
axi.text(0.04 * 2 * w1, 1.36 * hmax, "spectrum", color=GRAY, fontsize=6.5, ha="left", va="top")

# ================================================================ (c) one budget: the k-th sharpest rhythm
ax = axes_in(fig, 5.12, b0, 1.91, PH)
xs = np.logspace(np.log10(0.008), np.log10(0.5), 300)
YMIN, YMAX = 0.03, 40
ax.fill_between(xs, cap(xs), 2 * YMAX, color=LGRAY, lw=0, zorder=0)
# random networks as a density cloud (log-log): the regions holding 90% and 50% of their rhythms
H, xe_, ye_ = np.histogram2d(np.log10(Dc["rx"]), np.log10(Dc["rq"]), bins=150,
                             range=[[np.log10(0.01), np.log10(0.5)], [np.log10(YMIN), np.log10(YMAX)]])
H = gaussian_filter(H, 3.2)
srt = np.sort(H.ravel())[::-1]
cum = np.cumsum(srt) / srt.sum()
lev = [srt[np.searchsorted(cum, f)] for f in (0.9, 0.5)]
XC, YC = 10 ** (0.5 * (xe_[1:] + xe_[:-1])), 10 ** (0.5 * (ye_[1:] + ye_[:-1]))
ax.contourf(XC, YC, H.T, levels=[lev[0], lev[1], H.max() * 1.01], colors=CLOUD, zorder=1)
ax.contour(XC, YC, H.T, levels=[lev[0], lev[1]], colors=["#BDBDBD", "#A6A6A6"], linewidths=0.4, zorder=1.1)
ax.plot(xs, clk(xs), color=RING, lw=1.1, zorder=2)
ax.plot(Dc["chainx"], Dc["chainq"], color=NEST, lw=1.3, zorder=3, solid_capstyle="butt")
ax.plot(Dc["ringx"], Dc["ringq"], ls="none", marker="o", ms=2.7, mfc=ONE, mec="white", mew=0.4, zorder=4,
        clip_on=False)
ax.plot(xs, cap(xs), color=WEYL, lw=1.1, zorder=5)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.01, 0.5); ax.set_ylim(YMIN, YMAX)
ax.set_xticks([0.01, 0.1]); ax.set_xticklabels([r"$0.01$", r"$0.1$"])
ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels([r"$0.1$", r"$1$", r"$10$"])
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r"Rank $k/(d-1)$", labelpad=1.5)
ax.set_ylabel(r"Sharpness $Q_k$", labelpad=1.0)
along(ax, 0.22, cap(0.22) * 1.75, "capacity", (0.1, cap(0.1)), (0.4, cap(0.4)), color=WEYL)
ax.text(0.27, 20, "forbidden", color=FORBID_TEXT, fontsize=7, ha="center", va="center")
along(ax, 0.033, clk(0.033) / 1.62, "standalone clocks", (0.015, clk(0.015)), (0.045, clk(0.045)), color=RING)
ax.text(0.0115, 4.0, "one ring", color=ONE, fontsize=7, ha="left", va="bottom")
ax.text(0.0115, 1.0, "nested cycles", color=NEST, fontsize=7, ha="left", va="bottom")
ax.text(0.0115, 0.10, "random\nnetworks", color="#5C5C5C", fontsize=7, ha="left", va="bottom", linespacing=1.0)
xh = 0.040                                      # standalone clocks stay about a factor of two below the capacity
tip = 10 ** (0.6 / (PH * 72 / np.log10(YMAX / YMIN)))   # arrow tips end at the edges of the two lines
arrow(ax, (xh, clk(xh) * tip), (xh, cap(xh) / tip), "#4D4D4D", lw=0.6, ms=3.5)
xl = xh * 1.25                                  # its label inside the band, along the lines
along(ax, xl, 0.995 * np.sqrt(clk(xl) * cap(xl)), r"$\times2$", (xl / 1.3, cap(xl / 1.3)), (xl * 1.3, cap(xl * 1.3)),
      color="#4D4D4D")
header(fig, ax, "(c)", "One budget")

save(fig, os.path.join(HERE, "out", "fig2.pdf"))
print("fig2 size in", round(W, 3), round(Hh, 3))
