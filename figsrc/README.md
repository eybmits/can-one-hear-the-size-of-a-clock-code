# Figure sources

Python scripts that draw the four figures of the Letter and compute every data set they show, in the vocabulary
of the text: rhythms, sharpness Q, standalone clocks, one ring, nested cycles, random networks,
distinguishable and blurred rhythms, drive. The networks and parameters behind every panel are listed in the table at
the end of the Supplemental Material (Sec. VII, Code and data).

| File | Content |
|---|---|
| `style.py` | Shared style. Fonts are Computer Modern set by LaTeX, as in the REVTeX text. Figures are drawn at their printed width (column 246 pt, text 510 pt). Vivid, print-friendly palette with fixed meanings: blue (#1565C0) for standalone clocks and rhythms; rose (#D0287A) for one ring and its overtones; orange (#F07C12) for nested cycles; purple (#8E24AA) for drive; teal (#00A38C) for distinguishable rhythms; near-black (#1A1A1A) for bounds; gray for single-rhythm bounds, traffic, and forbidden regions; slate (#263845) for data. |
| `data.py` | Computes the data sets and stores them in `data/*.npz`, named by figure and panel: the ring of Fig. 2(b): return probability and occupation of its states at the times of the dials (`fig2b_ring`), the sharpness of the k-th sharpest rhythm in Fig. 2(c) (`fig2c_rank`: 440 random networks, a ring of 32 states, nested cycles of 1500 states), the density of Fig. 3(a) (`fig3a_density`), the bulk-chain curves of Fig. 3(b) (`fig3b_spacing`, 600 values of h per block size s), and the cost curves for Q > 10 and the simulated trace of Fig. 4 (`fig4a_costs`, `fig4b_trace`). Fig. 1 is a sketch and needs no data. |
| `fig1.py` | Fig. 1 (column width, page 1): three ways to hold many rhythms on ten states each: standalone clocks (one rhythm per ring), one ring (overtones of one cycle), and nested cycles (one rhythm per shortcut); darker states are slower. |
| `fig2.py` | Fig. 2 (text width, three panels): a drum (the ℓ-th mode of a square drum of side π against ℓ, Weyl's law ω_ℓ = (2/π)(1 + √(1 + πℓ)), and one mode as inset; the modes are computed in the script); a fading rhythm (the return probability of a ring of 20 states with its envelope, a window of dials whose hand, the mean of e^{2πij/d}, turns at ω and shrinks at a, and a window with the spectral line); one budget, the sharpness Q_k of the k-th sharpest rhythm against k/(d − 1) with the capacity of Eq. (2), on which standalone clocks (analytic, Q = cot(π/m) for rings of m states) and one ring stay a factor of two below and nested cycles come closest; random networks as a density cloud. |
| `fig3.py` | Fig. 3 (column width, two panels): Weyl law for clocks (standalone clocks; one large ring, arccot(Q)/π rhythms per state, exactly half the capacity for Q > 1; nested cycles); the price of a dip (s = 4 and 32). Both curves end at the same height, 10 % above the cost of standalone clocks: beyond it the damping angle for large s oscillates with h, so the cost stops being a smooth function of the spacing. |
| `fig4.py` | Fig. 4 (column width, two panels): states and drive per rhythm relative to standalone clocks; hearing a hidden clock. |
| `checkfig.py` | Overlap check: with `CHECK_OVERLAPS=1`, each script reports every label that touches other ink (lines, markers, other labels, tick marks); `CHECK_PAD` sets the margin in pixels at 300 dpi (default 2). |
| `bulkfast.py` | Bulk rhythms of nested cycles, found as zeros of the self-similar comparison function G (SM, Sec. II B). |
| `bulk.py` | Arbitrary-precision (mpmath) version of the same solver, used to validate `bulkfast.py`. |
| `data_fig4b.json` | Simulated cross-correlation, error bars and fit of the hidden machine in Fig. 4(b). |

## Fig. 4(a) units

States and drive per rhythm are divided by the cost of standalone clocks, N(n₀ − 1) + 1 states for N rhythms
(n₀ = 32 for Q > 10, about πQ states per rhythm), and by 2πQ of drive per rhythm, the least drive of one rhythm.
The capacity (sharp form) and the discrete form of the drive law both settle near one half. Of the single-rhythm
bounds, F_tot > 2πQ falls off like 1/N, while d ≥ max(n₀, 2N + 1) levels off at 2/(n₀ − 1) ≈ 0.065.

## Fig. 4(b) data

The points in `data_fig4b.json` were read losslessly from the vector file of the original figure; they are the
authors' simulated data. The three rhythms were fitted to the original fitted curves as damped cosines, with an rms
error of 1e-8. The full fit adds three real exponentials, rms error 4e-8. If the raw output of the inference protocol
is available (`data/worked_example.json` of the reproduction package), it can replace this file.

## Usage

```
python3 data.py            # about 10 min; most of it is the distinguishability curves of Fig. 3(b)
python3 fig1.py; python3 fig2.py; python3 fig3.py; python3 fig4.py
CHECK_OVERLAPS=1 python3 fig2.py   # optional: report labels that touch other ink
cp out/fig1.pdf out/fig2.pdf out/fig3.pdf out/fig4.pdf ../figs/
```

Requirements: Python 3 with numpy, scipy, mpmath, and matplotlib 3.8–3.10 (3.11.0 drops minus signs from LaTeX labels in PDFs), plus LaTeX with the `type1cm` and `cm-super`
packages for `text.usetex`.
