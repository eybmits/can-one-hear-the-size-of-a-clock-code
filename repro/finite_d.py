"""How fast is the sharp constant approached? Finite-size counts and the cycle-affinity price of theta chains.

SM Secs. II C and IV F. For Q0 in {2, 5, 10} and d up to 1500 (double-precision spectra of theta chains are reliable
there; see sharing_resolution.py accuracy), we count the pairs with Q > Q0 of the best theta chain (nested cycles) on a
parameter grid and compare with
  - independent clocks (standalone clocks): identical uniform rings glued at a hub (harmonics counted), best ring size;
  - the sharp form of the capacity (Corollary S1), printed as "Eq. (2)": N < 2(d-1) arccot(Q0)/pi.
For the best chain at each (Q0, d) we also give the spacing r = 2 Q tanh(sh/2) of its bulk rhythms, and the total cycle
affinity of its uniformly dressed version at eps = 3 (F_Sigma, and the total drive F~_Sigma = F_tot of SM Sec. IV B)
against that of independent rings with the same per-step affinity.

Usage:  python finite_d.py        (about 3 min); the output is finite_d.log
"""
import math
import sys

import numpy as np
from scipy.optimize import brentq

from models import theta_W, dress, stationary


def optimal_split(Js, T):
    """tau_c = sqrt(J_c^2 + 2 J_c/mu) with sum_c tau_c = T (total drive F_tot, SM Sec. IV B); a one-way edge
    (T = sum_c J_c) forces tau_c = J_c"""
    Js = np.asarray(Js, float)
    if len(Js) == 1:
        return np.array([max(T, Js[0])])
    if T <= Js.sum() * (1 + 1e-12):
        return Js.copy()
    f = lambda mu: np.sum(np.sqrt(Js ** 2 + 2 * Js / mu)) - T  # noqa: E731
    lo = 0.5 * 2 * np.sum(np.sqrt(Js)) ** 2 / T ** 2          # f(lo) > 0
    hi = 1.0
    while f(hi) > 0:
        hi *= 10
    mu = brentq(f, lo, hi, xtol=1e-300, rtol=1e-14, maxiter=500)
    return np.sqrt(Js ** 2 + 2 * Js / mu)


def qualities(W):
    ev = np.linalg.eigvals(W)
    up = ev[ev.imag > 1e-12 * np.abs(ev).max()]
    return up.imag / -up.real


def ring_count(d, Q0, t=1.0):
    """best glued identical uniform rings (per-step affinity eps with t = tanh(eps/2)); harmonics counted"""
    best = 0
    for n in range(3, d + 1):
        per = sum(1 for m in range(1, n) if 2 * m < n and t / math.tan(math.pi * m / n) > Q0)
        if per:
            best = max(best, ((d - 1) // (n - 1)) * per)
    return best


def best_chain(d, Q0):
    s0 = next(s for s in range(2, 1000) if 1 / math.tan(math.pi / (2 * s)) > Q0)
    best = (0, None)
    for s in range(s0, s0 + 9):
        B = d // s
        if B < 3:
            continue
        for sh in [0.01, 0.02, 0.04, 0.08, 0.15]:
            for a in [0.2, 0.5]:
                N = int((qualities(theta_W(s, sh / s, a, B)) > Q0).sum())
                if N > best[0]:
                    best = (N, (s, sh, a, B))
    return best


def affinity_costs(s, sh, a, B, eps):
    t = math.tanh(eps / 2)
    W = dress(theta_W(s, sh / s, a, B), t)
    d = len(W)
    p = stationary(W)
    F = W * p[None, :]
    np.fill_diagonal(F, 0)
    T, Jm = F + F.T, F - F.T
    rets = [s * (l + 1) - 1 for l in range(B)]
    Jl = np.array([Jm[0, r] for r in rets])
    FS = sum(eps * (r + 1) for r in rets)
    Ft = eps * B
    for j in range(d - 1):
        cs = [i for i, r in enumerate(rets) if r >= j + 1]
        Js = Jl[cs]
        ts = optimal_split(Js, T[j + 1, j])
        Ft += sum(math.log((tv + Jc) / (tv - Jc)) for tv, Jc in zip(ts, Js))
    return FS, Ft, qualities(W)


def main():
    for Q0 in [2.0, 5.0, 10.0]:
        print(f"== Q0 = {Q0}: bound density 2 arccot(Q0)/pi = {2 * math.atan(1 / Q0) / math.pi:.4f}", flush=True)
        for d in [100, 200, 400, 800, 1500]:
            N, par = best_chain(d, Q0)
            Nr = ring_count(d, Q0)
            Nb = 2 * (d - 1) * math.atan(1 / Q0) / math.pi
            line = f"   d={d:5d}: best theta chain N={N:4d}"
            if par:
                s, sh, a, B = par
                Qb = np.median(np.sort(qualities(theta_W(s, sh / s, a, B)))[::-1][:max(1, N)])
                r = 2 * Qb * math.tanh(sh / 2)
                line += f" (s={s}, sh={sh}, a={a}; median Q of counted pairs {Qb:.2f}, spacing r={r:.2f})"
            line += f"; independent rings N={Nr:4d}; ratio {N / max(Nr, 1):.3f}; Eq. (2) allows N < {Nb:.1f}"
            print(line, flush=True)
    print("== cycle affinity of dressed chains (eps = 3) vs independent rings at the same eps")
    eps = 3.0
    t = math.tanh(eps / 2)
    for (s, sh, a, B, Q0) in [(8, 0.24, 0.5, 25, 2.5), (8, 0.24, 0.5, 45, 2.5), (8, 0.24, 0.5, 90, 2.5)]:
        FS, Ft, Q = affinity_costs(s, sh, a, B, eps)
        N = int((Q > Q0).sum())
        n_ring = next(n for n in range(3, 10 ** 5) if t / math.tan(math.pi / n) > Q0)
        print(f"   s={s}, sh={sh}, B={B} (d={s * B}): N(Q>{Q0}) = {N}; F_Sigma = {FS:.0f}, F~_Sigma = {Ft:.0f}; independent "
              f"rings with eps = 3: {N} x {n_ring} states, F_Sigma = {N * n_ring * eps:.0f}; Eq. (7): F_Sigma > "
              f"{math.pi * N * (Q0 - 1 / Q0):.0f}", flush=True)


if __name__ == "__main__":
    main()
