"""Bulk modes of theta chains from the self-similar comparison function G (SM, Sec. II).

G(e^x) = (p/H) e^{kappa(x+H)} sum_j c_j Z^j,  Z = exp(-2 pi i x / H),
c_j = psihat(-2 pi j / H),
psihat(xi) = q^{kappa+i xi} * pi / sin(pi(kappa+i xi)) * P(-q^{1+kappa+i xi}) / P(-q),
P(y) = prod_{r>=0} (1 + y q^r) = (-y; q)_inf.
A bulk eigenvalue is lambda = e^x with sum_j c_j Z^j = 0; its damping angle is Im x - pi/2.
Neighbouring bulk modes have the same angle and |lambda| ratio e^H.
"""
import mpmath as mp

def coeffs(s, h, kappa, J, dps=60):
    mp.mp.dps = dps
    h = mp.mpf(h); kappa = mp.mpf(kappa)
    q = mp.e ** (-h); H = s * h
    Pq = mp.qp(q, q)  # (q;q)_inf = P(-q)
    c = {}
    for j in range(-J, J + 1):
        xi = -2 * mp.pi * j / H
        z = kappa + 1j * xi
        num = mp.qp(q ** (1 + z), q)   # P(-q^{1+z}) = (q^{1+z}; q)_inf
        c[j] = q ** z * mp.pi / mp.sin(mp.pi * z) * num / Pq
    return c, H

def bulk_mode(s, h, kappa=0.5, J=3, dps=60):
    c, H = coeffs(s, h, kappa, J, dps)
    F = lambda x: sum(c[j] * mp.exp(-2j * mp.pi * j * x / H) for j in c)
    # two-term start: c_0 + c_1 Z = 0  ->  x = i H/(2 pi) log(-c_0/c_1)
    x0 = 1j * H / (2 * mp.pi) * mp.log(-c[0] / c[1])
    x = mp.findroot(F, x0, tol=mp.mpf(10) ** (-(dps // 2)), maxsteps=200)
    if not (0 < x.imag < mp.pi):
        raise RuntimeError("root outside the upper half plane")
    delta = x.imag - mp.pi / 2
    return dict(delta=float(delta), Q=float(mp.cot(delta)), H=float(H), x=x,
                two_term_delta=float(x0.imag - mp.pi / 2))

if __name__ == "__main__":
    import time
    for s, h in ((4, 0.2), (4, 0.1), (4, 0.05), (4, 0.02)):
        t = time.time(); r = bulk_mode(s, h)
        print(s, h, "delta/(pi/2s)=%.4f Q=%.4f  2tanh(H/2)Q=%.4f (e^H-1)Q=%.4f (1-e^-H)Q=%.4f  y=%.4f  [%.1fs]" % (
            r["delta"] / (3.141592653589793 / (2 * s)), r["Q"], 2 * mp.tanh(r["H"] / 2) * r["Q"], (mp.e ** r["H"] - 1) * r["Q"], (1 - mp.e ** -r["H"]) * r["Q"],
            s / (3.141592653589793 / r["delta"] - 1), time.time() - t))
