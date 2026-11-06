"""Exact-arithmetic certification of generator spectra (stand-alone; Python standard library + mpmath).

Pipeline for a generator W with rational rates:
  1. exact stationary law (Gauss-Jordan over the rationals) and exact t = max_e |J_e| / T_e;
  2. exact characteristic polynomial det(sI - W) (Hessenberg reduction + recurrence over the rationals);
  3. roots to high precision (mpmath), rounded to rational centres z_i;
  4. inclusion discs: for monic p of degree n with distinct centres, w_i = p(z_i)/prod_{j!=i}(z_i - z_j) and every
     zero lies in the union of the discs D(z_i, n|w_i|); pairwise disjoint discs contain exactly one zero each
     (Lagrange interpolation and a homotopy w -> s*w). p(z_i), the products and |w_i|^2 are computed EXACTLY in
     Gaussian rationals; radii are rigorous rational upper bounds;
  5. for a disc in the open upper-left quadrant, Q = Im z/(-Re z) is enclosed in [(Im z-R)/(-Re z+R), (Im z+R)/(-Re z-R)].
"""
from fractions import Fraction as Fr
import json
import os

import mpmath as mp

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def load_rates(name):
    """returns d and the exact rate matrix R[dst][src] from a data file with 'edges' [src, dst] and decimal 'rates'."""
    D = json.load(open(os.path.join(DATA, name)))
    d = D["d"]
    R = [[Fr(0)] * d for _ in range(d)]
    for (src, dst), r in zip(D["edges"], D["rates"]):
        R[dst][src] += Fr(r)
    return d, R, D


def generator(R):
    d = len(R)
    W = [[R[i][j] if i != j else Fr(0) for j in range(d)] for i in range(d)]
    for j in range(d):
        W[j][j] = -sum(W[i][j] for i in range(d) if i != j)
    return W


def stationary(W):
    d = len(W)
    M = [W[i][:] + [Fr(0)] for i in range(d)]
    M[d - 1] = [Fr(1)] * d + [Fr(1)]
    for c in range(d):
        piv = next(r for r in range(c, d) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(d):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    p = [M[i][d] for i in range(d)]
    assert all(x > 0 for x in p), "stationary law not positive"
    assert all(sum(W[i][j] * p[j] for j in range(d)) == 0 for i in range(d))
    return p


def t_exact(W, p):
    """t = max_e |J_e|/T_e = tanh(eps/2) over edges with positive traffic (exact rational); None if some edge is one-way."""
    d = len(W)
    t = Fr(0)
    for i in range(d):
        for j in range(i + 1, d):
            f, b = W[i][j] * p[j], W[j][i] * p[i]
            if f + b > 0:
                if f == 0 or b == 0:
                    return None
                t = max(t, abs(f - b) / (f + b))
    return t


def charpoly(W):
    """exact coefficients of det(sI - W), highest degree first (Hessenberg reduction over the rationals)."""
    n = len(W)
    H = [row[:] for row in W]
    for m in range(1, n - 1):
        piv = next((i for i in range(m, n) if H[i][m - 1] != 0), None)
        if piv is None:
            continue
        if piv != m:
            H[m], H[piv] = H[piv], H[m]
            for row in H:
                row[m], row[piv] = row[piv], row[m]
        for i in range(m + 1, n):
            if H[i][m - 1] != 0:
                u = H[i][m - 1] / H[m][m - 1]
                H[i] = [a - u * b for a, b in zip(H[i], H[m])]
                for row in H:
                    row[m] += u * row[i]
    # recurrence for the characteristic polynomial of an upper Hessenberg matrix (coefficients low -> high)
    P = [[Fr(1)]]
    for m in range(1, n + 1):
        prev = P[m - 1]
        cur = [Fr(0)] + prev[:]                       # s * p_{m-1}
        for k, c in enumerate(prev):
            cur[k] -= H[m - 1][m - 1] * c              # - h_mm * p_{m-1}
        tprod = Fr(1)
        for i in range(1, m):
            tprod *= H[m - i][m - i - 1]
            coef = tprod * H[m - i - 1][m - 1]
            if coef != 0:
                for k, c in enumerate(P[m - i - 1]):
                    cur[k] -= coef * c
        P.append(cur)
    return P[n][::-1]


class GQ:
    __slots__ = ("re", "im")

    def __init__(self, re, im=Fr(0)):
        self.re, self.im = Fr(re), Fr(im)

    def __add__(self, o): return GQ(self.re + o.re, self.im + o.im)
    def __sub__(self, o): return GQ(self.re - o.re, self.im - o.im)
    def __mul__(self, o): return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    def abs2(self): return self.re * self.re + self.im * self.im


def sqrt_upper(x, digits=40):
    with mp.workdps(digits + 10):
        s = Fr(mp.nstr(mp.sqrt(mp.mpf(x.numerator) / x.denominator) * (1 + mp.mpf(10) ** (-digits)), digits + 5))
    while s * s < x:
        s *= Fr(1000001, 1000000)
    return s


def certify_roots(c, dps=130, cdigits=110):
    """c: exact monic coefficients, highest first. Returns [(centre GQ, radius Fraction)], disjointness flag."""
    n = len(c) - 1
    with mp.workdps(dps):
        roots = mp.polyroots([mp.mpf(x.numerator) / x.denominator for x in c], maxsteps=4000, extraprec=4 * dps)
        Z = [GQ(Fr(mp.nstr(mp.re(z), cdigits)), Fr(mp.nstr(mp.im(z), cdigits))) for z in roots]

    def peval(z):
        acc = GQ(0)
        for coef in c:
            acc = acc * z + GQ(coef)
        return acc
    out = []
    for i, z in enumerate(Z):
        num, den = peval(z), GQ(1)
        for j, y in enumerate(Z):
            if j != i:
                den = den * (z - y)
        w2 = num.abs2() / den.abs2()
        out.append((z, n * sqrt_upper(w2) if w2 > 0 else Fr(0)))
    disjoint = all((out[i][0] - out[j][0]).abs2() > (out[i][1] + out[j][1]) ** 2
                   for i in range(n) for j in range(i))
    return out, disjoint


def certified_qualities(W):
    """rigorous lower bounds on the qualities of all non-real pairs (descending), and the largest disc radius."""
    c = charpoly(W)
    assert c[-1] == 0, "zero is not an eigenvalue"
    c = c[:-1]                                         # divide by s
    discs, disjoint = certify_roots(c)
    assert disjoint, "inclusion discs overlap"
    lows = []
    for z, R in discs:
        if z.im > R:
            assert -z.re > R
            lows.append((z.im - R) / (-z.re + R))
    return sorted(lows, reverse=True), max(float(R) for _, R in discs)


def ivlog_ratio(a, b):
    """interval enclosure of ln(a/b) for positive rationals a, b."""
    iv = mp.iv
    return iv.log(iv.mpf(a.numerator) / a.denominator) - iv.log(iv.mpf(b.numerator) / b.denominator)
