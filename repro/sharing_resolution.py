"""The price of sharing: state saving of theta chains versus the spectral resolution of their rhythms.

Fig. 3(b) and SM Sec. III B. A theta chain T(s,h,a,B), the chain of nested cycles T(s,h,kappa,B) of SM Sec. II with
a = kappa, carries one rhythm per block of s states. Deep in the bulk its spectrum is log-periodic: the rhythm
of the next block is the present one rescaled by e^{-sh}. Two consecutive rhythms lambda and lambda e^{-sh}, with quality Q,
are therefore separated in frequency by r = 2 Q tanh(sh/2) times their mean half-width (a dip between two equal Lorentzians
needs r > 2/sqrt(3)). The bulk rhythm is a zero of the renewal function Lambda of Eq. (S12), located here by Newton's method
in the middle of a chain with B*s*h >= 40, so that edge effects are negligible.

The state cost per rhythm is s for the chain, pi/arccot(Q) - 1 for standalone clocks of the same quality (rings glued at a
hub; continuum form, a lower bound on the integer ring cost), and pi/(2 arccot Q) for the bound of Corollary S1.

Part 2 (bounded rate spread): at a fixed span e^{h(d-1)} of all rates, the number of pairs with Q > Q0 of explicit chains
grows in proportion to d, because the rhythms crowd more densely in frequency.

Usage:  python sharing_resolution.py        (about 7 min); writes ../data/sharing_resolution.json
        python sharing_resolution.py span   (part 2 only, a few seconds)
        python sharing_resolution.py accuracy   (double-precision eigenvalues vs Newton-refined zeros of Lambda)
"""
import json
import math
import os
import sys

import numpy as np

from models import theta_W


def renewal(lam, s, h, a, B):
    """Lambda(lam) and Lambda'(lam) for T(s,h,a,B), evaluated with log-sums (Eq. (S12))"""
    p = 1 - math.exp(-a * s * h)
    j = np.arange(s * B)
    x = lam * np.exp(h * j)                     # lam q^{-j}
    lg = np.cumsum(np.log1p(x))
    dl = np.cumsum(np.exp(h * j) / (1 + x))
    al = p * (1 - p) ** np.arange(B)
    al[-1] = (1 - p) ** (B - 1)
    idx = s * (np.arange(B) + 1) - 1
    P = np.exp(-lg[idx])
    return 1 - np.sum(al * P), np.sum(al * P * dl[idx])


def bulk_rhythm(s, sh, a):
    """a zero of Lambda on the l = 0 branch in the middle of a long chain"""
    h = sh / s
    B = max(60, int(math.ceil(40 / sh)))
    target = 2 * math.pi * h * math.exp(-sh * B / 2)       # 2 pi / (traversal time of the cycle through block B/2)
    best = None
    for m in [0.3, 1.0, 3.0]:
        for dlt in np.linspace(0.6, 2.2, 9) * math.pi / (2 * s):
            z = m * target * np.exp(1j * (math.pi / 2 + dlt))
            for _ in range(200):
                L, dL = renewal(z, s, h, a, B)
                step = L / dL
                if abs(step) > 0.3 * abs(z):
                    step *= 0.3 * abs(z) / abs(step)
                z -= step
                if abs(step) < 1e-14 * abs(z):
                    break
            ang = np.angle(z) - math.pi / 2
            if abs(renewal(z, s, h, a, B)[0]) < 1e-9 and 0 < ang < math.pi / s and abs(math.log(abs(z) / target)) < 3:
                if best is None or ang < best[1]:
                    best = (z, ang)
    if best is None:
        return None
    Q = 1 / math.tan(best[1])
    return dict(s=s, sh=sh, a=a, B=B, Q=Q, r=2 * Q * math.tanh(sh / 2))


def check_against_eig(s, sh, a, B=60):
    """the bulk zero agrees with the median bulk eigenvalue of an explicit chain (consistency check)"""
    W = theta_W(s, sh / s, a, B)
    ev = np.linalg.eigvals(W)
    up = ev[ev.imag > 1e-12 * abs(ev).max()]
    Q = np.sort(up.imag / (-up.real))[::-1][:B]
    return float(np.median(Q[B // 4: 3 * B // 4]))


def bounded_span(Q0=2.0, s=5, a=0.5, spans=(2.0, 4.0)):
    """pairs with Q > Q0 of T(s,h,a,B) at fixed rate span e^{h(d-1)}, as d grows"""
    print(f"bounded rate span: s={s}, a={a}, Q0={Q0}; densities N/(d-1) against 1/s = {1 / s:.3f}, "
          f"rings {1 / (next(n for n in range(3, 99) if 1 / math.tan(math.pi / n) > Q0) - 1):.3f}, "
          f"bound {2 / math.pi * math.atan(1 / Q0):.3f}")
    for lnR in spans:
        for d in [100, 200, 400, 800]:
            B = d // s
            h = lnR / (s * B - 1)
            ev = np.linalg.eigvals(theta_W(s, h, a, B))
            up = ev[ev.imag > 1e-12 * abs(ev).max()]
            N = int((up.imag / (-up.real) > Q0).sum())
            print(f"   span e^{lnR:.0f} = {math.exp(lnR):6.1f}, d = {s * B:4d}: N = {N:4d}, density {N / (s * B - 1):.3f}")


def eig_accuracy():
    """double-precision eigenvalues of explicit chains vs Newton-refined zeros of the renewal function"""
    worst = 0.0
    for s, sh, a, B in [(4, 0.1, 0.35, 120), (8, 0.01, 0.35, 160), (8, 0.05, 0.35, 160), (8, 0.3, 0.35, 80),
                        (16, 0.02, 0.35, 80), (16, 0.1, 0.35, 120), (4, 0.04, 0.5, 375)]:
        h = sh / s
        ev = np.linalg.eigvals(theta_W(s, h, a, B))
        up = ev[ev.imag > 1e-12 * abs(ev).max()]
        up = up[np.angle(up) - math.pi / 2 < math.pi / s]          # the coherent branch used in the analysis
        for z0 in up:
            z = z0
            for _ in range(60):
                L, dL = renewal(z, s, h, a, B)
                step = L / dL
                z -= step
                if abs(step) < 1e-15 * abs(z):
                    break
            worst = max(worst, abs(z - z0) / abs(z0))
        print(f"   s={s:2d} sh={sh} B={B} (d={s * B}): {len(up)} pairs on the branch delta < pi/s refined", flush=True)
    print(f"largest relative shift of a double-precision eigenvalue under Newton refinement on Lambda: {worst:.1e}")


def main():
    ring_cost = lambda Q: math.pi / math.atan(1 / Q) - 1      # noqa: E731
    bound_cost = lambda Q: math.pi / (2 * math.atan(1 / Q))  # noqa: E731
    out = []
    for s in [4, 8, 16, 32]:
        for sh in np.geomspace(0.01, 1.5, 40):
            res = None
            for a in [0.2, 0.35, 0.5]:
                r = bulk_rhythm(s, float(sh), a)
                if r is not None and (res is None or r["Q"] > res["Q"]):
                    res = r
            if res is None:
                continue
            res["cost_ratio"] = s / ring_cost(res["Q"])            # chain states per rhythm / ring states per rhythm
            res["bound_ratio"] = bound_cost(res["Q"]) / ring_cost(res["Q"])
            res["ring_int_cost"] = next(n for n in range(3, 10 ** 6) if 1 / math.tan(math.pi / n) > res["Q"]) - 1
            out.append(res)
            print(f"s={s:2d} sh={sh:.4f} a={res['a']:.2f}  Q={res['Q']:7.3f}  r={res['r']:6.3f}  "
                  f"chain/ring states = {res['cost_ratio']:.3f} (integer rings: {s / res['ring_int_cost']:.3f}; "
                  f"bound/ring = {res['bound_ratio']:.3f})", flush=True)
    for s, sh, a in [(4, 0.1, 0.35), (8, 0.1, 0.35), (8, 0.3, 0.35), (16, 0.1, 0.35)]:
        z = bulk_rhythm(s, sh, a)
        B = int(math.ceil(12 / sh))
        print(f"check s={s} sh={sh}: bulk zero Q = {z['Q']:.4f}; median bulk eigenvalue Q of an explicit chain with "
              f"B={B} blocks (d={s * B}) = {check_against_eig(s, sh, a, B):.4f}")
    for s in [4, 8, 16, 32]:
        rows = [o for o in out if o["s"] == s]
        print(f"s={s}: bulk rhythm on the l = 0 branch found up to r = {max(o['r'] for o in rows):.3f} "
              f"(sh = {max(o['sh'] for o in rows):.3f}); beyond, it leaves the branch (Q < cot(pi/s))")
        res = [o for o in rows if o["r"] >= 2 / math.sqrt(3)]
        if res:
            print(f"s={s}: among resolved rhythms (r >= 2/sqrt3) the smallest chain/ring state ratio is "
                  f"{min(o['cost_ratio'] for o in res):.3f}")
        print(f"s={s}: smallest chain/ring state ratio overall {min(o['cost_ratio'] for o in rows):.3f} "
              f"at r = {min(rows, key=lambda o: o['cost_ratio'])['r']:.3f}")
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "sharing_resolution.json")
    json.dump(out, open(path, "w"), indent=1)


if __name__ == "__main__":
    np.seterr(all="ignore")          # Newton steps from poor starting points may diverge; such candidates are rejected
    if sys.argv[1:] == ["accuracy"]:
        eig_accuracy()
        sys.exit()
    if sys.argv[1:] != ["span"]:
        main()
    bounded_span()
