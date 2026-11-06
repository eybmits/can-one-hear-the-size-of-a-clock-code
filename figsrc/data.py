"""Compute (or extract) every data set shown in the figures and store it in data/*.npz.

Run:  python3 data.py            (all)   or   python3 data.py fig2 fig3 fig4
"""
import json, os, sys
import numpy as np
from scipy.linalg import eig
from scipy.optimize import minimize_scalar, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)


# --------------------------------------------------------------------------- generators
def theta_chain(s, h, kappa, B):
    """T(s,h,kappa,B): d=sB states, exit rate q^j, return to the hub (state 0) w.p. p at block ends."""
    d = s * B
    q = np.exp(-h)
    p = 1.0 - np.exp(-kappa * s * h)
    W = np.zeros((d, d))
    for j in range(d):
        r = q ** j
        if j == d - 1:
            W[0, j] += r
        elif (j + 1) % s == 0:
            W[0, j] += p * r
            W[j + 1, j] += (1 - p) * r
        else:
            W[j + 1, j] += r
        W[j, j] -= r
    return W


def ring(n, f=1.0, b=0.0):
    W = np.zeros((n, n))
    for j in range(n):
        W[(j + 1) % n, j] += f
        W[(j - 1) % n, j] += b
        W[j, j] -= f + b
    return W


def pairs(W, tol=1e-9):
    """Upper-half-plane eigenvalues (one per conjugate pair)."""
    ev = np.linalg.eigvals(W)
    scale = np.abs(ev).max()
    return ev[ev.imag > tol * scale]


def damping_angles(W):
    lam = pairs(W)
    return np.sort(np.arctan2(-lam.real, lam.imag))  # delta = arctan(a/omega), increasing


def stationary(W):
    w, v = eig(W)
    k = np.argmin(np.abs(w))
    p = np.real(v[:, k]); p = p / p.sum()
    return p


# --------------------------------------------------------------------------- Figure 2 (budget)
def fig2_ring_return(d=20, rate=1.0, dial_t=(0, 10, 20, 40)):
    """Many copies of a one-way ring of d states, started together in state 0 (Fig. 2b):
    return probability P_00(t), and the occupation of the states at the times of the dials (left inset),
    p_j(t) = (1/d) sum_m exp(lambda_m t) exp(-2 pi i m j / d). The right inset is the Lorentzian line of the sharpest
    rhythm, drawn from its frequency w1 and decay rate a1."""
    m = np.arange(d)
    lam = rate * (np.exp(2j * np.pi * m / d) - 1.0)
    t = np.linspace(0, 80, 4001)
    P = np.real(np.exp(np.outer(t, lam)).sum(axis=1)) / d
    a1, w1 = -lam[1].real, lam[1].imag
    dial_t = np.array(dial_t, dtype=float)
    dial_p = np.real(np.exp(np.outer(dial_t, lam)) @ np.exp(-2j * np.pi * np.outer(m, m) / d)) / d
    assert np.allclose(dial_p.sum(axis=1), 1) and dial_p.min() > -1e-12
    np.savez(os.path.join(OUT, "fig2b_ring.npz"), t=t, P=P, d=d, a1=a1, w1=w1, dial_t=dial_t,
             dial_p=np.clip(dial_p, 0, None))
    print(f"  fig2b ring d={d}: omega={w1:.4f}, a={a1:.4f}, Q={w1/a1:.3f}; dials at t={[int(x) for x in dial_t]}")


def random_generator(d, rng, density=0.6, spread=2.0):
    W = np.zeros((d, d))
    mask = rng.random((d, d)) < density
    np.fill_diagonal(mask, False)
    W[mask] = np.exp(spread * rng.standard_normal(mask.sum()))
    # make it strongly connected by adding a directed cycle through a random permutation
    perm = rng.permutation(d)
    for a, b in zip(perm, np.roll(perm, -1)):
        W[b, a] += np.exp(spread * rng.standard_normal())
    np.fill_diagonal(W, 0.0)
    W -= np.diag(W.sum(axis=0))
    return W


def fig2_rank():
    """Sharpness Q_k of the k-th sharpest rhythm against k/(d-1) (Fig. 2c): 40 random networks for each odd d from
    5 to 25, a one-way ring of 32 states, and the chain T(4, 0.01, 1/2, 375) of 1500 states, the best chain of
    SM Sec. II C at Q0 = 2 (344 rhythms with Q > 2). The capacity, Eq. (2), is Q_k < 2(d-1)/(pi k)."""
    rng = np.random.default_rng(11)
    rx, rq = [], []
    for d in range(5, 26, 2):
        for _ in range(40):
            W = random_generator(d, rng, density=0.5, spread=2.0)
            dl = damping_angles(W)
            if len(dl) == 0:
                continue
            rx.append(np.arange(1, len(dl) + 1) / (d - 1)); rq.append(1 / np.tan(dl))
    rx, rq = np.concatenate(rx), np.concatenate(rq)
    n = 32
    k = np.arange(1, n // 2)
    ringx, ringq = k / (n - 1), 1 / np.tan(np.pi * k / n)
    dC = damping_angles(theta_chain(4, 0.01, 0.5, 375))
    chainx, chainq = np.arange(1, len(dC) + 1) / 1499, 1 / np.tan(dC)
    np.savez(os.path.join(OUT, "fig2c_rank.npz"), rx=rx, rq=rq, ringx=ringx, ringq=ringq, chainx=chainx, chainq=chainq)
    frac = lambda x, q: np.max(q * x) / (2 / np.pi)  # noqa: E731  share of the capacity reached
    print(f"  fig2c random rhythms {len(rx)} (best {frac(rx, rq):.0%} of the capacity), ring {frac(ringx, ringq):.0%}, "
          f"chain {frac(chainx, chainq):.0%} with {np.sum(chainq > 2)} rhythms Q > 2")


# --------------------------------------------------------------------------- Figure 3 (sharing)
def fig3_density():
    Q0 = np.exp(np.linspace(np.log(0.3), np.log(40), 4000))
    bound = Q0 * np.minimum(0.5, 2 / np.pi * np.arctan(1 / Q0))
    chain = np.empty_like(Q0); rings = np.empty_like(Q0)
    for i, q in enumerate(Q0):
        s = 2
        while 1 / np.tan(np.pi / (2 * s)) <= q:
            s += 1
        chain[i] = q / s
        n0 = 3
        while 1 / np.tan(np.pi / n0) <= q:
            n0 += 1
        rings[i] = q / (n0 - 1)
    np.savez(os.path.join(OUT, "fig3a_density.npz"), Q0=Q0, bound=bound, chain=chain, rings=rings)


def fig3_resolution():
    from bulkfast import bulk_mode
    out = {}
    for s in (4, 8, 16, 32):
        def best(h):
            r = minimize_scalar(lambda k: bulk_mode(s, h, k)["delta"], bounds=(0.03, 0.97),
                                method="bounded", options={"xatol": 1e-6})
            return r.fun, r.x
        # choose h so that the spacing runs from ~0.05 to where the cost exceeds ~1.3; sampled densely, because for
        # large s the damping angle starts to oscillate with h near the end of the range (see fig3.py)
        hs = np.exp(np.linspace(np.log(0.04 * np.pi / (2 * s * s)), np.log(10.0 * np.pi / (2 * s * s)), 600))
        X, Y, D, K, H = [], [], [], [], []
        for h in hs:
            dlt, kap = best(h)
            Q = 1 / np.tan(dlt); Hs = s * h
            x = 2 * np.tanh(Hs / 2) * Q
            y = s / (np.pi / dlt - 1)
            X.append(x); Y.append(y); D.append(dlt); K.append(kap); H.append(h)
            if y > 1.30 or x > 1.6:
                break
        X, Y, D = np.array(X), np.array(Y), np.array(D)
        out[f"s{s}_x"] = X; out[f"s{s}_y"] = Y; out[f"s{s}_delta"] = D
        out[f"s{s}_kappa"] = np.array(K); out[f"s{s}_h"] = np.array(H)
        out[f"s{s}_bound"] = 1 / (2 - 2 * D / np.pi)
        i = int(np.argmax(Y >= 1.0))
        xc = X[i - 1] + (1.0 - Y[i - 1]) * (X[i] - X[i - 1]) / (Y[i] - Y[i - 1])
        mono = int(np.argmin(np.r_[True, (np.diff(X) > 0) & (np.diff(Y) > 0)])) or len(X)
        print(f"  fig3b s={s}: x {X[0]:.3f}..{X[-1]:.3f}, cost crosses rings at spacing {xc:.3f}, "
              f"rises monotonically up to x={X[mono - 1]:.3f} (y={Y[mono - 1]:.3f}), kappa* {K[0]:.2f}..{K[-1]:.2f}")
    np.savez(os.path.join(OUT, "fig3b_spacing.npz"), **out)


# --------------------------------------------------------------------------- Figure 4 (costs)
def fig4_costs():
    Q0 = 10.0
    ac = np.arctan(1 / Q0)
    N = np.arange(1, 31)
    def dmin_states(n, t):
        d = n + 1
        while True:
            rhs = np.sum(np.arctan(t / np.tan(np.pi * (2 * np.arange(1, n + 1) - 1) / (2 * (d - 1)))))
            if n * np.arctan(Q0) < rhs:
                return d
            d += 1
    out = {"N": N}
    for lab, t in (("inf", 1.0), ("one", np.tanh(0.5))):
        multi = np.array([dmin_states(n, t) for n in N])
        n0 = 3
        while t / np.tan(np.pi / n0) <= Q0:
            n0 += 1
        out["multi_" + lab] = multi / N
        out["single_" + lab] = np.maximum(n0, 2 * N + 1) / N     # sharpest rhythm and the count of rhythms
        out["rings_" + lab] = (N * (n0 - 1) + 1) / N
        out["n0_" + lab] = n0
    # affinity: discrete form  N arctan Q0 < sum_{m<=N} arctan(F/(2 pi m))
    Fmin = []
    for n in N:
        g = lambda F: np.sum(np.arctan(F / (2 * np.pi * np.arange(1, n + 1)))) - n * np.arctan(Q0)
        Fmin.append(brentq(g, 1e-6, 1e6))
    Fmin = np.array(Fmin)
    out["aff_multi"] = Fmin / N
    out["aff_single"] = 2 * np.pi * Q0 / N
    out["aff_rings"] = np.full(N.shape, 2 * np.pi * Q0)
    np.savez(os.path.join(OUT, "fig4a_costs.npz"), **out)
    print(f"  fig4a N=5: d>= {out['multi_inf'][4]*5:.0f} (eps=inf), {out['multi_one'][4]*5:.0f} (eps=1); "
          f"single {out['single_inf'][0]:.0f}/{out['single_one'][0]:.0f}; rings {out['rings_inf'][4]*5:.0f}/{out['rings_one'][4]*5:.0f}; "
          f"F>= {Fmin[4]:.2f}")


def fig4_trace():
    src = json.load(open(os.path.join(HERE, "data_fig4b.json")))
    np.savez(os.path.join(OUT, "fig4b_trace.npz"), tau=np.array(src["tau"]), C=np.array(src["C"]),
             lo=np.array(src["lo"]), hi=np.array(src["hi"]),
             comps=np.array([[c["A"], c["a"], c["w"], c["ph"]] for c in src["components"]]),
             exps=np.array(src["exps"]))


if __name__ == "__main__":
    # data sets by figure of the Letter; Fig. 1 is a sketch and needs no data
    which = sys.argv[1:] or ["fig2", "fig3", "fig4"]
    if "fig2" in which:            # Fig. 2: fading ring, sharpness of the k-th rhythm
        print("fig2"); fig2_ring_return(); fig2_rank()
    if "fig3" in which:            # Fig. 3: Weyl law for clocks, the price of a dip
        print("fig3"); fig3_density(); fig3_resolution()
    if "fig4" in which:            # Fig. 4: states and drive per rhythm, the simulated trace
        print("fig4"); fig4_costs(); fig4_trace()
