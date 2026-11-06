"""Shared figure style for the Letter.

* Fonts: Computer Modern through LaTeX, as in the REVTeX text (cmr10).
* Sizes: figures are drawn at their printed width (column 246 pt, text 510 pt), so nothing is rescaled
  by \\includegraphics; axis labels 8 pt, tick labels and annotations 7 pt, panel labels 8.5 pt.
* Colors: one vivid, print-friendly palette (blue/rose/orange/purple/teal, safe for common color-vision
  deficiencies), used with the same meaning in every figure:
    standalone clocks, rhythms           -> RING (blue)
    one ring and its overtones           -> ONE (rose)
    nested cycles                        -> NEST (orange)
    drive                                -> TURN (purple)
    distinguishable rhythms              -> DIST (teal)
    bounds and totals                    -> INK (charcoal);  single-rhythm bounds, traffic -> GRAY
    forbidden regions                    -> light gray fill;  data -> DATA (slate)
"""
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

PT = 1 / 72.27                      # TeX point in inches
COLW = 246.0 * PT                   # \columnwidth
TEXTW = 510.0 * PT                  # \textwidth

# Vivid, print-friendly palette (Material-like hues; safe for common color-vision deficiencies)
RING = "#1565C0"          # standalone clocks, rhythms of probabilities
RING_DARK = "#0D3C7A"
RING_LIGHT = "#5FA2EE"
RING_FILL = "#D3E5FB"
ONE = "#D0287A"           # one ring and its overtones
ONE_DARK = "#8F1450"
ONE_LIGHT = "#F08BB5"
ONE_FILL = "#FBE1EE"
NEST = "#F07C12"          # nested cycles
NEST_DARK = "#B14E00"
NEST_LIGHT = "#FFBF7A"
NEST_FILL = "#FFE4C8"
TURN = "#8E24AA"          # circulation, roots of the spanning-tree polynomial, drive
TURN_FILL = "#F0DDF6"
DIST = "#00A38C"          # distinguishable rhythms
DIST_TEXT = "#00796B"
DIST_FILL = "#D0F1EB"
INK = "#1A1A1A"           # bounds and totals
WEYL = "#C62828"          # Weyl's law (Fig. 2a) and its counterpart for clocks, the capacity (Fig. 2c)
DATA = "#263845"          # data
GRAY = "#707070"          # single-rhythm bounds, traffic
MGRAY = "#ADADAD"
LGRAY = "#ECECEC"         # forbidden regions
FORBID_TEXT = "#8C8C8C"
BLACK = INK
SEQ_RING = [RING_DARK, RING, RING_LIGHT]                    # three rhythms of the hidden machine

RC = {
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{bm}",
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.labelpad": 2.0,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 2.8,
    "ytick.major.size": 2.8,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.pad": 2.5,
    "ytick.major.pad": 2.0,
    "lines.linewidth": 1.2,
    "lines.markersize": 3.0,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "legend.handletextpad": 0.5,
    "legend.borderaxespad": 0.3,
    "legend.labelspacing": 0.25,
    "hatch.linewidth": 0.5,
    "pdf.fonttype": 42,
    "savefig.dpi": 600,
    "path.simplify": False,
}
plt.rcParams.update(RC)


def along(ax, x, y, text, p, q, **kw):
    """Text at (x, y), rotated to the screen direction of the data segment p -> q."""
    import numpy as np
    P = ax.transData.transform(np.array([p, q], dtype=float))
    ang = np.degrees(np.arctan2(P[1, 1] - P[0, 1], P[1, 0] - P[0, 0]))
    kw.setdefault("rotation_mode", "anchor")
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    kw.setdefault("fontsize", 7)
    return ax.text(x, y, text, rotation=ang, **kw)


def figure(width, height):
    return plt.figure(figsize=(width, height))


def axes_in(fig, x0, y0, w, h, **kw):
    """Axes at (x0, y0) with size (w, h), all in inches from the lower-left corner."""
    W, H = fig.get_size_inches()
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H], **kw)


def panel_label(fig, ax, label, dx=0.0, dy=0.0, x=None):
    """Panel label '(a)' just above the upper-left corner of the axes frame.
    x (inches) fixes the horizontal position, so labels of one column line up."""
    W, H = fig.get_size_inches()
    bb = ax.get_position()
    xx = (x / W) if x is not None else bb.x0 + dx / W
    fig.text(xx, bb.y1 + (dy + 0.03) / H, label, ha="left", va="bottom", fontsize=9, color=BLACK)


def minor(ax, x=True, y=True, n=2):
    if x:
        ax.xaxis.set_minor_locator(AutoMinorLocator(n))
    if y:
        ax.yaxis.set_minor_locator(AutoMinorLocator(n))


def forbidden(ax, x, y_lo, y_hi, text=None, tx=None, ty=None, **kw):
    ax.fill_between(x, y_lo, y_hi, color=LGRAY, lw=0, zorder=0)
    if text:
        ax.text(tx, ty, text, color=FORBID_TEXT, fontsize=7, ha=kw.get("ha", "left"), va=kw.get("va", "bottom"),
                transform=kw.get("transform", ax.transData))


def save(fig, path):
    import os
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if os.environ.get("CHECK_OVERLAPS"):
        from checkfig import report_overlaps
        report_overlaps(fig, os.path.basename(path), pad_px=int(os.environ.get("CHECK_PAD", "2")))
    fig.savefig(path, metadata={"Creator": None, "Producer": None, "CreationDate": None, "Title": None})
    plt.close(fig)


def header(fig, ax, label, title=None, dy=0.04, margin=0.02):
    """Panel label and title, '(a) Title', centered above the panel's axes box (its original box, so that
    equal-aspect drawings keep the row's baseline). The text is shifted inward if it would leave the figure."""
    W, H = fig.get_size_inches()
    bb = ax.get_position(original=True)
    xc, y0 = (bb.x0 + bb.x1) / 2 * W, bb.y1 * H + dy
    text = label + (r"\enspace " + title if title else "")
    t = fig.text(xc / W, y0 / H, text, ha="center", va="baseline", fontsize=8, color=INK)
    r = fig._get_renderer()
    ext = t.get_window_extent(renderer=r)
    lo, hi = margin * fig.dpi, (W - margin) * fig.dpi
    shift = (hi - ext.x1) if ext.x1 > hi else ((lo - ext.x0) if ext.x0 < lo else 0.0)
    if shift:
        t.set_position((xc / W + shift / (W * fig.dpi), y0 / H))
    return t
