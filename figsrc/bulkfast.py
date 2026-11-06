"""Fast double-precision bulk-mode solver for theta chains (see bulk.py for the formulas).

log c_j = z log q + log(pi) - log sin(pi z) + log P(-q^{1+z}) - log P(-q),  z = kappa - 2 pi i j / H,
with q^{i xi} = exp(2 pi i j / s) a root of unity, so
log P(-q^{1+z}) = sum_r log(1 - e^{-h(1+kappa+r)} e^{2 pi i j/s}).
The root of sum_j c_j Z^j is found in the normalised variable W = (c_1/c_0) Z near W = -1.
"""
import numpy as np

def _logP(h, kappa, j, s):
    nr = int(60.0 / h) + 10
    r = np.arange(nr, dtype=float)
    a = np.exp(-h * (1.0 + kappa + r)) * np.exp(2j * np.pi * j / s)
    return np.sum(np.log1p(-a))

def _logsinpi(z):
    # log sin(pi z) for z = kappa + i y, stable for large |y|
    y = z.imag
    if y > 0:
        # sin(pi z) = -e^{-i pi z} (1 - e^{2 i pi z}) / (2i),  e^{2 i pi z} small
        return -1j * np.pi * z + np.log(-(1 - np.exp(2j * np.pi * z)) / (2j))
    else:
        # sin(pi z) = e^{i pi z} (1 - e^{-2 i pi z}) / (2i),  e^{-2 i pi z} small
        return 1j * np.pi * z + np.log((1 - np.exp(-2j * np.pi * z)) / (2j))

def logc(s, h, kappa, J=3):
    H = s * h
    lq = -h
    logPq = np.sum(np.log1p(-np.exp(-h * (1.0 + np.arange(int(60.0 / h) + 10)))))
    out = {}
    for j in range(-J, J + 1):
        z = kappa - 2j * np.pi * j / H
        out[j] = z * lq + np.log(np.pi) - _logsinpi(z) + _logP(h, kappa, j, s) - logPq
    return out, H

def bulk_mode(s, h, kappa=0.5, J=3):
    lc, H = logc(s, h, kappa, J)
    # r_j = c_j c_0^{j-1} / c_1^j
    lr = {j: lc[j] + (j - 1) * lc[0] - j * lc[1] for j in lc}
    def F(W):
        return sum(np.exp(lr[j]) * W ** j for j in lr if np.real(lr[j]) > -700)
    def dF(W):
        return sum(j * np.exp(lr[j]) * W ** (j - 1) for j in lr if np.real(lr[j]) > -700)
    W = -1.0 + 0j
    for _ in range(100):
        step = F(W) / dF(W); W -= step
        if abs(step) < 1e-15: break
    logZ = np.log(W) + lc[0] - lc[1]
    imx = H / (2 * np.pi) * np.real(logZ)
    delta = imx - np.pi / 2
    return dict(delta=delta, Q=1 / np.tan(delta), H=H, W=W)

if __name__ == "__main__":
    for s, h in ((4, 0.2), (4, 0.1), (4, 0.05), (4, 0.02), (8, 0.05), (8, 0.02)):
        r = bulk_mode(s, h)
        print(s, h, "delta/(pi/2s)=%.6f Q=%.6f |W+1|=%.2e" % (r["delta"] / (np.pi / (2 * s)), r["Q"], abs(r["W"] + 1)))
