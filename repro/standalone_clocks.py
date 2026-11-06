"""Standalone clocks: N one-way rings of n0 states that share one state (SM Sec. I).

Ring c (c = 1 slowest, ..., N fastest) has n0 - 1 own states, each left at rate v_c; the shared state enters ring c at
rate nu_c. With equal entry rates (nu_c = v_c) only the fastest ring keeps its sharpness; with the graded rates
    nu_N = v_N,   nu_c = v_c [1 + (n0 - 1) sum_{c' > c} nu_{c'} / v_{c'}],
every ring tends to Q = cot(pi/n0) as the speed ratios grow. The script prints the damping angle
delta = arctan(a/omega) of the sharpest rhythm near each ring's frequency, in 60-digit arithmetic.
"""
import mpmath as mp

mp.mp.dps = 60


def generator(n0, v, nu):
    N = len(v)
    d = 1 + N * (n0 - 1)
    W = mp.zeros(d, d)
    idx = lambda c, j: 1 + c * (n0 - 1) + (j - 1)          # own state j = 1..n0-1 of ring c
    for c in range(N):
        W[idx(c, 1), 0] += nu[c]                          # shared state -> first own state
        for j in range(1, n0 - 1):
            W[idx(c, j + 1), idx(c, j)] += v[c]
        W[0, idx(c, n0 - 1)] += v[c]                      # last own state -> shared state
    for k in range(d):
        W[k, k] = -sum(W[i, k] for i in range(d) if i != k)
    return W


def graded(n0, v):
    N = len(v)
    nu = [mp.mpf(0)] * N
    nu[N - 1] = v[N - 1]
    for c in range(N - 2, -1, -1):
        nu[c] = v[c] * (1 + (n0 - 1) * sum(nu[k] / v[k] for k in range(c + 1, N)))
    return nu


def ring_angles(n0, v, nu):
    ev = mp.eig(generator(n0, v, nu), left=False, right=False)
    out = []
    for c, vc in enumerate(v):
        # rhythms of ring c live near |lambda| ~ v_c; pick the sharpest one in that decade
        cand = [l for l in ev if mp.im(l) > 0 and vc / 30 < abs(l) < 30 * vc]
        out.append(min(mp.atan(-mp.re(l) / mp.im(l)) for l in cand))
    return out


if __name__ == "__main__":
    n0, N = 8, 3
    print("pi/n0 =", mp.nstr(mp.pi / n0, 10))
    for ratio in (mp.mpf(10) ** 2, mp.mpf(10) ** 4, mp.mpf(10) ** 6, mp.mpf(10) ** 8):
        v = [ratio ** c for c in range(N)]                # slowest first
        eq = ring_angles(n0, v, v)
        gr = ring_angles(n0, v, graded(n0, v))
        print(f"ratio 1e{int(mp.log10(ratio))}:  equal entry {[mp.nstr(x, 5) for x in eq]}   "
              f"graded entry {[mp.nstr(x, 7) for x in gr]}   max |graded - pi/n0| = "
              f"{mp.nstr(max(abs(x - mp.pi / n0) for x in gr), 3)}")
