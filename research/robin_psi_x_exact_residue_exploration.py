"""Exploratory exact-residue study of psi_X=f_X'/F_X for the kappa=4 Robin benchmark.

This is a high-precision numerical exploration, not an interval certificate.
The residues are evaluated from the exact half/full characteristic functions,
not from a finite-product approximation.
"""
import json
import math
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq

mp.mp.dps = 60
KAPPA = mp.mpf(4)
W1, W2, WR = mp.mpf(25) / 66, mp.mpf(4) / 33, mp.mpf(1) / 2
D1, D2, DR = mp.mpf(25), mp.mpf(1) / 4, mp.mpf(1)
N = 60


def mat_mp(k, h, D):
    z = k * h / mp.sqrt(D)
    c, s = mp.cos(z), mp.sin(z)
    return ((c, s / (mp.sqrt(D) * k)),
            (-mp.sqrt(D) * k * s, c))


def matmul(A, B):
    return tuple(tuple(sum(A[i][q] * B[q][j] for q in range(2))
                       for j in range(2)) for i in range(2))


def left_state(k):
    M = matmul(mat_mp(k, W2, D2), mat_mp(k, W1, D1))
    return (M[0][0] + M[0][1] * KAPPA,
            M[1][0] + M[1][1] * KAPPA)


def full_characteristic(k):
    y, p = left_state(k)
    M = mat_mp(k, WR, DR)
    yr = M[0][0] * y + M[0][1] * p
    pr = M[1][0] * y + M[1][1] * p
    return pr + KAPPA * yr


def right_characteristic(k):
    return mp.cos(k / 2) + KAPPA * mp.sin(k / 2) / k


def mat_float(k, h, D):
    z = k * h / math.sqrt(D)
    c, s = math.cos(z), math.sin(z)
    return np.array([[c, s / (math.sqrt(D) * k)],
                     [-math.sqrt(D) * k * s, c]])


def left_state_float(k):
    M = mat_float(k, float(W2), float(D2)) @ mat_float(k, float(W1), float(D1))
    return M @ np.array([1.0, 4.0])


def full_characteristic_float(k):
    y, p = left_state_float(k)
    yr, pr = mat_float(k, 0.5, 1.0) @ np.array([y, p])
    return pr + 4.0 * yr


def right_characteristic_float(k):
    return math.cos(k / 2) + 4.0 * math.sin(k / 2) / k


def scan_roots(fun, kmax, step, count):
    grid = np.arange(1.0e-5, kmax + step, step)
    vals = np.array([fun(float(x)) for x in grid])
    roots = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if fa * fb < 0:
            x = brentq(fun, a, b, xtol=1.0e-13, rtol=1.0e-14)
            if not roots or x - roots[-1] > 1.0e-7:
                roots.append(x)
                if len(roots) == count:
                    break
    if len(roots) != count:
        raise RuntimeError(f"root scan found {len(roots)} roots, expected {count}")
    return roots


# Obtain root seeds, then refine each root at 60 decimal digits.
k_global_seed = scan_roots(full_characteristic_float, 400.0, 0.05, N + 2)
k_right_seed = [brentq(right_characteristic_float,
                       (2 * j - 1) * math.pi + 1.0e-9,
                       2 * j * math.pi - 1.0e-9,
                       xtol=1.0e-13, rtol=1.0e-14)
                for j in range(1, N + 2)]
k_global = [mp.findroot(full_characteristic,
                        (mp.mpf(x) - mp.mpf("1e-5"), mp.mpf(x) + mp.mpf("1e-5")))
             for x in k_global_seed]
k_right = [mp.findroot(right_characteristic,
                       (mp.mpf(x) - mp.mpf("1e-5"), mp.mpf(x) + mp.mpf("1e-5")))
            for x in k_right_seed]
lam = [x * x for x in k_global]
R = [x * x for x in k_right]
c = lam[1]
lam1, R1 = lam[0], R[0]

# P_G(z)=D_G(sqrt(-z))/D_G(0), P_R(z)=Y_R(sqrt(-z))/Y_R(0).
# For a root k, P_G'(-k^2)=-D_G'(k)/(2 k D_G(0)); D_G(0)=24.
def PG_prime_at(k):
    return -mp.diff(full_characteristic, k) / (2 * k * 24)


def PR_at(k):
    return right_characteristic(k) / 3


PG_prime_c = PG_prime_at(k_global[1])
PR_at_c = PR_at(mp.sqrt(c))
normalization = ((1 - c / R1) / PR_at_c
                 * (c * PG_prime_c) / (1 - c / lam1))

eta, rho, residues = [], [], []
for i in range(2, N + 2):  # half-spectrum index i=2,...,N+1
    lam_n = lam[i]        # full eigenvalue lambda_{i+1}
    eta_i = lam_n - c
    rho_i = R[i - 1] - c
    PG_prime_n = PG_prime_at(k_global[i])
    PR_n = PR_at(k_global[i])
    B_i = (normalization * PR_n / (1 - lam_n / R1)
           * (1 - lam_n / lam1) * (1 - lam_n / c) / PG_prime_n)
    eta.append(eta_i)
    rho.append(rho_i)
    residues.append(B_i)


def evaluate(t, count=N):
    t = mp.mpf(t)
    exps = [mp.exp(-eta[i] * t) for i in range(count)]
    f = mp.fsum(residues[i] * exps[i] for i in range(count))
    fp = mp.fsum(-eta[i] * residues[i] * exps[i] for i in range(count))
    fpp = mp.fsum(eta[i] ** 2 * residues[i] * exps[i] for i in range(count))
    F = 1 - mp.fsum(residues[i] / eta[i] * exps[i] for i in range(count))
    psi = fp / F
    Npsi = fpp * F - fp * f
    psi_prime = Npsi / F ** 2
    return {"f": f, "F": F, "f_prime": fp, "f_second": fpp,
            "psi": psi, "Npsi": Npsi, "psi_prime": psi_prime}


root = mp.findroot(lambda t: evaluate(t)["Npsi"],
                   (mp.mpf("0.014"), mp.mpf("0.015")))

sample_times = ["0.0015", "0.003", "0.005", "0.0075", "0.010",
                "0.012", "0.014", "0.0144", "0.0145", "0.015",
                "0.020", "0.035", "0.070", "0.100", "0.200", "0.500"]
samples = {t: {k: mp.nstr(evaluate(t)[k], 20)
                for k in ("f", "F", "Npsi", "psi_prime")}
           for t in sample_times}
truncation = {
    t: {str(n): mp.nstr(evaluate(t, n)["psi_prime"], 18)
        for n in (10, 15, 20, 30, 40, 50, 60)}
    for t in ("0.010", "0.0144", "0.0145", "0.020", "0.035", "0.070")
}

result = {
    "status": "exploratory_not_interval_certified",
    "parameters": {"w1": "25/66", "D1": "25", "w2": "4/33",
                   "D2": "1/4", "right_width": "1/2",
                   "D_right": "1", "kappa": "4"},
    "root_data": {
        "lambda1": mp.nstr(lam[0], 30),
        "lambda2_c": mp.nstr(c, 30),
        "lambda3": mp.nstr(lam[2], 30),
        "R1": mp.nstr(R1, 30),
        "R2": mp.nstr(R[1], 30),
        "eta2": mp.nstr(eta[0], 30),
        "rho2": mp.nstr(rho[0], 30),
        "B2": mp.nstr(residues[0], 30),
        "B3": mp.nstr(residues[1], 30)
    },
    "psi_prime_zero_candidate": mp.nstr(root, 40),
    "samples": samples,
    "truncation_stability": truncation,
    "runtime": {"python": platform.python_version(),
                "mpmath": mp.__version__, "numpy": np.__version__,
                "scipy": scipy.__version__, "platform": platform.platform(),
                "machine": platform.machine(), "cpu_count": __import__("os").cpu_count()}
}

result["truncation_stability_early"] = {
    t: {str(n): mp.nstr(evaluate(t, n)["psi_prime"], 18)
        for n in (10, 15, 20, 30, 40, 50, 60)}
    for t in ("0.0015", "0.003", "0.005", "0.0075")
}
result["root_truncation_stability"] = {
    str(n): mp.nstr(mp.findroot(lambda t: evaluate(t, n)["Npsi"],
                                (mp.mpf("0.014"), mp.mpf("0.015"))), 28)
    for n in (10, 15, 20, 30, 40, 50, 60)
}
print(json.dumps(result, indent=2))
