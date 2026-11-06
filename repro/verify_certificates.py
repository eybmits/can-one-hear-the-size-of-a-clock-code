"""Exact verification of the computer-assisted statements of the paper (SM Sec. IV E cites items (i) and (ii)).

  (i)   single driven cycle: the exactly stalled 27-state network (one driven 20-edge cycle plus a current-free side
        path) has 2*pi*3*Q_(3)/F_Sigma > 1.39, i.e. the additive affinity law fails;
  (i')  the original, unstalled 27-state theta graph: 2*pi*3*Q_(3)/F_Sigma > 1.098;
  (ii)  62-state comb: Q_(30) > T_30(62), i.e. the ring-mediant (per-index ring) conjecture fails;
  (iii) 11-state cactus: Q_(2) exceeds the second-largest per-cycle value, i.e. the per-index affinity law fails;
  (iv)  near-equilibrium networks W6, W8 (eps = 0.004): Q_(k)/t exceeds the uniform-ring value cot(pi k/d).

All printed lower bounds are truncated and all printed upper bounds are rounded up, so every printed inequality is
implied by the exact computation. Usage:  python verify_certificates.py        (about 1 min)
"""
import json
import math
import os
from fractions import Fraction as Fr

import mpmath as mp

from exactcert import (DATA, load_rates, generator, stationary, t_exact, certified_qualities, ivlog_ratio)

iv = mp.iv
iv.dps = 60


def _raw_to_fr(raw):
    sgn, man, exp, _ = raw
    v = Fr(int(man)) * Fr(2) ** exp
    return -v if sgn else v


def iv_lo(x):
    """exact rational lower endpoint of an mpmath interval"""
    return _raw_to_fr(x._mpi_[0])


def iv_hi(x):
    """exact rational upper endpoint of an mpmath interval"""
    return _raw_to_fr(x._mpi_[1])


def lo(F, k):
    """decimal string <= F with k digits after the point"""
    n = math.floor(F * 10 ** k)
    return f"{Fr(n, 10 ** k):.{k}f}" if n >= 0 else "-" + f"{Fr(-n, 10 ** k):.{k}f}"


def hi(F, k):
    """decimal string >= F with k digits after the point"""
    n = math.ceil(F * 10 ** k)
    return f"{Fr(n, 10 ** k):.{k}f}" if n >= 0 else "-" + f"{Fr(-n, 10 ** k):.{k}f}"


def path_affinity(R, P):
    return sum((ivlog_ratio(R[P[a + 1]][P[a]], R[P[a]][P[a + 1]]) for a in range(len(P) - 1)), iv.mpf(0))


def theta_paths(d):
    v = d - 1
    return [0] + list(range(1, 16)) + [v], [0, 16, 17, 18, v], [0] + list(range(19, 26)) + [v]


def additive_law_stalled():
    d, R, _ = load_rates("additive_law_witness_stalled_d27.json")
    W = generator(R)
    p = stationary(W)
    P16, P4, P8 = theta_paths(d)
    cur = lambda a, b: R[b][a] * p[a] - R[a][b] * p[b]  # noqa: E731
    assert all(cur(P8[i], P8[i + 1]) == 0 for i in range(len(P8) - 1)), "side path not stalled"
    J16, J4 = cur(P16[0], P16[1]), cur(P4[0], P4[1])
    assert J16 > 0 and J4 == -J16
    F = path_affinity(R, P16) - path_affinity(R, P4)       # cycle 0 -> P16 -> v -> P4 (reversed) -> 0
    lows, rmax = certified_qualities(W)
    Q3 = lows[2]
    ratio_lo = iv_lo(2 * iv.pi * 3 * iv.mpf(Q3.numerator) / Q3.denominator / F)
    eff_lo = Q3 * 3 / iv_hi(F)
    print(f"(i)   stalled d={d}: side-path current = 0 exactly; single current cycle of 20 edges, "
          f"F in [{lo(iv_lo(F), 13)}, {hi(iv_hi(F), 13)}]; Q_(1..3) >= {[lo(x, 6) for x in lows[:3]]}; "
          f"2 pi 3 Q_(3)/F >= {lo(ratio_lo, 4)}; Q0*N/F_Sigma >= {lo(eff_lo, 4)} as Q0 -> Q_(3)   "
          f"[additive law fails: {ratio_lo > 1}; max disc radius {rmax:.1e}]")


def additive_law():
    d, R, _ = load_rates("additive_law_witness_d27.json")
    W = generator(R)
    p = stationary(W)
    paths = list(theta_paths(d))
    J = [R[P[1]][P[0]] * p[P[0]] - R[P[0]][P[1]] * p[P[1]] for P in paths]
    A = [path_affinity(R, P) for P in paths]
    sgn = [1 if x > 0 else -1 for x in J]
    odd = [i for i in range(3) if sgn.count(sgn[i]) == 1][0]
    # conformal decomposition of a theta graph: the odd path closes a cycle with each of the other two paths
    F_sigma = sum(abs(A[odd] - A[i]) for i in range(3) if i != odd)
    lows, rmax = certified_qualities(W)
    Q3 = lows[2]
    ratio_lo = iv_lo(2 * iv.pi * 3 * iv.mpf(Q3.numerator) / Q3.denominator / F_sigma)
    eff = Fr("0.1501") * 3 / iv_hi(F_sigma)
    print(f"(i')  unstalled d={d}: path currents {[f'{float(x):+.3e}' for x in J]}; F_Sigma in "
          f"[{lo(iv_lo(F_sigma), 13)}, {hi(iv_hi(F_sigma), 13)}]; Q_(1..3) >= {[lo(x, 6) for x in lows[:3]]}; "
          f"2 pi 3 Q_(3)/F_Sigma >= {lo(ratio_lo, 4)}; Q0*N/F_Sigma >= {lo(eff, 4)} at Q0 = 0.1501   "
          f"[max disc radius {rmax:.1e}]")


def ring_mediant():
    d, R, D = load_rates("c4_counterexample_d62.json")
    k = D["k"]
    W = generator(R)
    lows, rmax = certified_qualities(W)
    num = sum((iv.sin(2 * iv.pi * m / d) for m in range(1, k + 1)), iv.mpf(0))
    den = sum((1 - iv.cos(2 * iv.pi * m / d) for m in range(1, k + 1)), iv.mpf(0))
    T_hi = iv_hi(num / den)
    Qk = lows[k - 1]
    print(f"(ii)  d={d}: Q_({k}) >= {lo(Qk, 14)} > {hi(T_hi, 14)} >= T_{k}({d}): {Qk > T_hi}   "
          f"[max disc radius {rmax:.1e}]")


def per_index_affinity():
    d, R, _ = load_rates("per_index_affinity_witness_d11.json")
    eps = Fr(1, 1000)
    for i in range(d):            # make every edge reversible: reverse rate = 1e-3 * forward rate
        for j in range(d):
            if i != j and R[i][j] > 0 and R[j][i] == 0:
                R[j][i] = eps * R[i][j]
    W = generator(R)
    lows, rmax = certified_qualities(W)
    cycles = [[0, 2, 3, 4, 5], [3, 6, 7], [6, 8, 9, 10]]     # the current-carrying cycles of the cactus
    vals = []
    for c in cycles:
        n = len(c)
        Fc = sum((ivlog_ratio(R[c[(a + 1) % n]][c[a]], R[c[a]][c[(a + 1) % n]]) for a in range(n)), iv.mpf(0))
        th = (iv.exp(Fc / n) - 1) / (iv.exp(Fc / n) + 1)          # tanh(F_c / 2 n_c)
        for m in range(1, n):
            if 2 * m < n:
                vals.append(iv_hi(th / iv.tan(iv.pi * m / n)))
    vals.sort(reverse=True)
    print(f"(iii) d={d}: Q_(2) >= {lo(lows[1], 6)} > {hi(vals[1], 6)} >= second-largest "
          f"tanh(F_c/2n_c)cot(pi m/n_c): {lows[1] > vals[1]}   [max disc radius {rmax:.1e}]")


def near_equilibrium():
    D = json.load(open(os.path.join(DATA, "near_equilibrium_witnesses.json")))
    for name, k in [("W6", 2), ("W8", 3)]:
        w = D[name]
        d = w["d"]
        R = [[Fr(0)] * d for _ in range(d)]
        for (src, dst), r in zip(w["edges"], w["rates"]):
            R[dst][src] += Fr(r)
        W = generator(R)
        p = stationary(W)
        t = t_exact(W, p)
        lows, rmax = certified_qualities(W)
        Qk = lows[k - 1]
        cot_hi = iv_hi(1 / iv.tan(iv.pi * k / d))            # t cot(pi k/d) <= t * cot_hi
        eps = 2 * mp.atanh(mp.mpf(t.numerator) / t.denominator)
        print(f"(iv)  {name}: d={d}, exact t = {lo(t, 15)}... (largest edge affinity {float(eps):.10f}); "
              f"Q_(1..{k}) >= {[lo(x, 14) for x in lows[:k]]}; Q_({k})/t >= {lo(Qk / t, 4)}; "
              f"Q_({k})/[t cot(pi {k}/{d})] >= {lo(Qk / (t * cot_hi), 4)}   [max disc radius {rmax:.1e}]")


if __name__ == "__main__":
    additive_law_stalled()
    additive_law()
    per_index_affinity()
    near_equilibrium()
    ring_mediant()
