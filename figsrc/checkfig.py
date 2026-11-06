"""Overlap check for the figures: does any text label touch other ink?

Run a figure script with CHECK_OVERLAPS=1; style.save() then calls `report_overlaps(fig)`.
For every visible text (annotations, legend entries, panel headers) the
figure is rendered once without that text. The pixels the text adds are dilated by `pad_px`, and
any ink of the remaining figure inside that region is reported. Light fills (white, the gray
forbidden regions, the coral band, the light-blue wedge) count as background, not ink.
"""
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.text import Text
from scipy.ndimage import binary_dilation

BACKGROUND = np.array([
    [255, 255, 255],      # paper
    [236, 236, 236],      # forbidden regions (LGRAY)
    [251, 233, 228],      # coral band, alpha 0.13 on white
    [238, 244, 250],      # light-blue wedge
    [227, 232, 243],      # proof wedge
], float)


def _render(fig, dpi):
    canvas = FigureCanvasAgg(fig)
    fig.set_dpi(dpi)
    canvas.draw()
    return np.asarray(canvas.buffer_rgba())[..., :3].astype(float)


def _texts(fig):
    out = []
    for t in fig.texts:
        out.append((t, "figure text"))
    for ax in fig.axes:
        for t in ax.texts:
            out.append((t, "annotation"))
        legends = [a for a in ax.get_children() if a.__class__.__name__ == "Legend"]
        for leg in legends:
            for t in leg.get_texts():
                out.append((t, "legend"))
        # tick and axis labels are not toggled (hiding one moves the axis label), but they count as ink
    return [(t, kind) for t, kind in out if isinstance(t, Text) and t.get_visible() and t.get_text().strip()]


def report_overlaps(fig, name="", dpi=300, pad_px=2, min_px=3):
    full = _render(fig, dpi)
    found = []
    for t, kind in _texts(fig):
        t.set_visible(False)
        without = _render(fig, dpi)
        t.set_visible(True)
        own = np.abs(full - without).max(axis=2) > 25
        if not own.any():
            continue
        if t.get_bbox_patch() is not None:          # a boxed label hides what lies under it
            continue
        region = binary_dilation(own, iterations=pad_px)
        dist = np.min(np.linalg.norm(without[..., None, :] - BACKGROUND[None, None], axis=-1), axis=-1)
        ink = dist > 40
        hits = region & ink
        n = int(hits.sum())
        if n >= min_px:
            ys, xs = np.nonzero(hits)
            found.append((n, kind, t.get_text(), xs.mean() / dpi, (full.shape[0] - ys.mean()) / dpi))
    label = f"[{name}] " if name else ""
    if not found:
        print(f"{label}no text overlaps")
    for n, kind, txt, x, y in sorted(found, key=lambda f: -f[0]):
        print(f"{label}{kind:12s} {txt!r:40s} touches {n:4d} px of ink at ({x:.2f} in, {y:.2f} in)")
    return found
