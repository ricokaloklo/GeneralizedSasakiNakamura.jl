# Independent reference: exact-extremal Kerr (a = M = 1) QNMs from the Richartz (PRD 93, 064062, 2016)
# five-term recurrence around r = 2M, with the corrected P3 = -2[... - 2im(1 + 6iω)].
# Convergent-subspace condition: backward recursion of two solutions with QR at every step, then the
# 2x2 boundary determinant of rows (3.2), (3.3). Angular eigenvalue by the spin-weighted spherical
# spectral method. Reproduces Tables I, V, VI of the paper (s = -2) to the printed six decimals.
# Usage: from richartz_extremal_reference import solve; solve(-2, l, m, seed, N=8000) -> (omega, A, |det|)
import numpy as np


def angular_A(s, m, c, lmax=40, target=None):
    lmin = max(abs(m), abs(s))
    ls = np.arange(lmin, lmax + 3)
    n = len(ls)
    X = np.zeros((n, n))
    for i, l in enumerate(ls):
        X[i, i] = -m * s / (l * (l + 1)) if l > 0 else 0.0
        if i + 1 < n:
            L = l + 1
            v = np.sqrt((L * L - m * m) * (L * L - s * s)) / (L * np.sqrt((2 * l + 1) * (2 * l + 3)))
            X[i, i + 1] = X[i + 1, i] = v
    X2 = (X @ X)
    k = lmax - lmin + 1
    H = np.diag([l * (l + 1) - s * (s + 1) for l in ls[:k]]).astype(complex) - c * c * X2[:k, :k] + 2 * c * s * X[:k, :k]
    ev = np.linalg.eigvals(H)
    if target is None:
        return ev
    return ev[np.argmin(abs(ev - target))]




def five_coeffs(s, m, w, lam, n):
    I = 1j
    P1 = 4 * (I * m - I * w - s)
    P2 = 2 * (1 - 4 * I * m + 16 * I * w)
    P3 = -2 * (1 + 2 * s + 2 * lam + 2 * I * w * (4 + 23 * I * w) - 2 * I * m * (1 + 6 * I * w))
    P4 = 4 * (I * m - I * w + s)
    P5 = -4 * (I * m - I * w + s) * (1 + 4 * I * w)
    P6 = -3 - 8 * I * w
    P7 = 2 * (1 + 6 * I * w - 8 * w * w)
    return n * n + n, P1 * n, -2 * n * n + P2 * n + P3, P4 * n + P5, n * n + P6 * n + P7


def boundary_det(s, m, w, lam, N):
    n = np.arange(N + 2, dtype=float)
    al, be, ga, de, ep = five_coeffs(s, m, w, lam, n)
    # window holds (a_{k+1}, a_k, a_{k-1}, a_{k-2}) for two solutions, as columns
    W = np.array([[1.0, 0.0], [0.0, 1.0], [0.3, -0.2], [-0.1, 0.4]], dtype=complex)
    for k in range(N, 2, -1):
        new = -(al[k] * W[0] + be[k] * W[1] + ga[k] * W[2] + de[k] * W[3]) / ep[k]
        W = np.vstack([W[1:], new])
        q, _ = np.linalg.qr(W)
        W = q
    # W rows = (a3, a2, a1, a0)
    a3, a2, a1, a0 = W
    E1 = al[1] * a2 + be[1] * a1 + ga[1] * a0
    E2 = al[2] * a3 + be[2] * a2 + ga[2] * a1 + de[2] * a0
    M = np.array([E1, E2])
    return np.linalg.det(M) / (np.linalg.norm(M[0]) * np.linalg.norm(M[1]))


def solve(s, l, m, w0, lam0=None, N=4000, tol=1e-13, it=60):
    lp = [lam0 if lam0 is not None else l * (l + 1) - s * (s + 1)]
    def F(w):
        lam = angular_A(s, m, w, lmax=max(40, l + 30), target=lp[0]); lp[0] = lam
        return boundary_det(s, m, w, lam, N)
    # secant iteration
    x0, x1 = complex(w0), complex(w0) + 1e-4
    f0, f1 = F(x0), F(x1)
    for _ in range(it):
        if f1 == f0: break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) > 0.02: x2 = x1 + 0.02 * (x2 - x1) / abs(x2 - x1)
        x0, f0, x1 = x1, f1, x2
        f1 = F(x1)
        if abs(x1 - x0) < tol: break
    return x1, lp[0], abs(f1)
