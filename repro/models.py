"""Network models used in the paper (stand-alone; numpy only).

Convention: W[i, j] is the rate j -> i, columns sum to zero.
  theta_W(s, h, a, B)                 theta chain T(s,h,a,B) = nested cycles of SM Sec. II, a = kappa (all edges one-way)
  theta_return_paths(s, h, a, B, r)   the same with every return edge replaced by a fast return path of r states
  dress(W, t)                         keep stationary law and traffic, multiply every current by t (inverse undressing)
  kinesin_W(row, ATP, ADP, P, F)      six-state kinesin network of Liepelt and Lipowsky, PRL 98, 258102 (2007)
"""
import math

import numpy as np


def theta_W(s, h, a, B):
    d = s * B
    q = math.exp(-h)
    p = 1 - math.exp(-a * s * h)
    W = np.zeros((d, d))
    for j in range(d):
        mu = q ** j
        if j == d - 1:
            W[0, j] += mu
        elif j % s == s - 1:
            W[0, j] += p * mu
            W[j + 1, j] += (1 - p) * mu
        else:
            W[j + 1, j] += mu
        W[j, j] -= mu
    return W


def theta_return_paths(s, h, a, B, r, R_fast=50.0):
    if r == 0:
        return theta_W(s, h, a, B)
    q = math.exp(-h)
    p = 1 - math.exp(-a * s * h)
    n_main = s * B
    d = n_main + r * B
    W = np.zeros((d, d))
    nxt = n_main
    for j in range(n_main):
        mu = q ** j
        if j == n_main - 1:
            ret, fwd = mu, 0.0
        elif j % s == s - 1:
            ret, fwd = p * mu, (1 - p) * mu
        else:
            ret, fwd = 0.0, mu
        if fwd > 0:
            W[j + 1, j] += fwd
        if ret > 0:
            path = list(range(nxt, nxt + r))
            nxt += r
            W[path[0], j] += ret
            for u, v in zip(path[:-1], path[1:]):
                W[v, u] += R_fast
            W[0, path[-1]] += R_fast
    W -= np.diag(W.sum(0))
    return W


def stationary(W):
    ev, vec = np.linalg.eig(W)
    p = np.real(vec[:, np.argmin(abs(ev))])
    return p / p.sum()


def dress(W, t):
    p = stationary(W)
    F = W * p[None, :]
    np.fill_diagonal(F, 0)
    T, J = F + F.T, F - F.T
    Ft = (T + t * J) / 2
    np.fill_diagonal(Ft, 0)
    Wt = Ft / p[None, :]
    return Wt - np.diag(Wt.sum(0))


def qualities(W):
    ev = np.linalg.eigvals(W)
    up = ev[ev.imag > 1e-12 * np.max(np.abs(ev))]
    return np.sort(up.imag / (-up.real))[::-1]


# ------------------------------------------------------------------------------------------------- kinesin
KEQ = 4.9e11  # uM, balance condition of Liepelt and Lipowsky
# Table I rows (k25, k52, k12 [1/(uM s)], k21, k56 = k61, k65 printed, k16 [1/(uM s)]); force parameters theta, chi12,
# chi56 from the appendix (published for the Carter-Cross and Visscher sets only; irrelevant at F = 0)
ROWS = {
    "CarterCross": dict(k25=3e5, k52=0.24, k12=2.0, k21=100, k56=100, k65p=0.02, k16=0.02,
                        theta=0.65, chi12=0.25, chi56=0.15),
    "Visscher": dict(k25=3e5, k52=0.24, k12=1.8, k21=100, k56=200, k65p=0.09, k16=0.02,
                     theta=0.3, chi12=0.25, chi56=0.05),
    "Schief": dict(k25=3e5, k52=0.24, k12=4.0, k21=100, k56=225, k65p=0.25, k16=0.02,
                   theta=0.65, chi12=0.25, chi56=0.15),
}


def kinesin_W(row, ATP, ADP, P, F_pN=0.0):
    r = ROWS[row]
    k65 = r["k25"] * r["k12"] * r["k56"] ** 2 / (r["k52"] * r["k21"] * r["k16"] * KEQ)
    k54 = (r["k52"] / r["k25"]) ** 2 * r["k21"]
    Fb = F_pN / 0.5                                   # l F / kBT with kBT/l = 0.5 pN
    phi_c12 = 2 / (1 + math.exp(r["chi12"] * Fb))
    phi_c56 = 2 / (1 + math.exp(r["chi56"] * Fb))
    th = r["theta"]
    s = lambda i: i - 1                               # states 1..6 -> 0..5
    rates = {
        (s(1), s(2)): r["k12"] * ATP * phi_c12, (s(2), s(1)): r["k21"] * phi_c12,
        (s(2), s(5)): r["k25"] * math.exp(-th * Fb), (s(5), s(2)): r["k52"] * math.exp((1 - th) * Fb),
        (s(5), s(6)): r["k56"] * phi_c56, (s(6), s(5)): k65 * ADP * phi_c56,
        (s(6), s(1)): r["k56"] * phi_c56, (s(1), s(6)): r["k16"] * P * phi_c56,
        (s(4), s(5)): r["k12"] * ATP * phi_c12, (s(5), s(4)): k54 * phi_c12,
        (s(2), s(3)): r["k56"] * phi_c56, (s(3), s(2)): k65 * ADP * phi_c56,
        (s(3), s(4)): r["k56"] * phi_c56, (s(4), s(3)): r["k16"] * P * phi_c56,
    }
    W = np.zeros((6, 6))
    for (i, j), rate in rates.items():
        W[j, i] += rate
    return W - np.diag(W.sum(0))
