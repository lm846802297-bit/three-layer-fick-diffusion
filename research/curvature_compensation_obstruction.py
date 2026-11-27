"""Interval checks for a Robin spectral family approaching L1=R1=lambda_2.

Only mpmath interval arithmetic is required.
"""
import json
import platform
import subprocess
from fractions import Fraction
from pathlib import Path
from mpmath import iv

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "curvature_compensation_obstruction_input.json").read_text())
iv.dps = DATA["interval_precision_decimal_digits"]
I = iv.mpf


class Dual:
    """First-order interval dual number."""
    def __init__(self, value, derivative=0):
        self.v = I(value)
        self.d = I(derivative)

    def __add__(self, other):
        other = as_dual(other)
        return Dual(self.v + other.v, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-as_dual(other))

    def __rsub__(self, other):
        return as_dual(other) + (-self)

    def __mul__(self, other):
        other = as_dual(other)
        return Dual(self.v * other.v, self.d * other.v + self.v * other.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_dual(other)
        return Dual(self.v / other.v,
                    (self.d * other.v - self.v * other.d) / (other.v * other.v))

    def __rtruediv__(self, other):
        return as_dual(other) / self


def as_dual(value):
    return value if isinstance(value, Dual) else Dual(value)


def dsin(x):
    return Dual(iv.sin(x.v), iv.cos(x.v) * x.d)


def dcos(x):
    return Dual(iv.cos(x.v), -iv.sin(x.v) * x.d)


def dtan(x):
    return Dual(iv.tan(x.v), x.d / (iv.cos(x.v) ** 2))


def crossing_derivative(k_bounds):
    k = Dual(I(k_bounds), 1)
    dw1 = Dual(I(25) / 66)
    dw2 = Dual(I(4) / 33)
    d1 = Dual(5)
    d2 = Dual(I(1) / 2)
    z1 = k * dw1 / d1
    z2 = k * dw2 / d2
    m1 = [
        [dcos(z1), dsin(z1) / (d1 * k)],
        [-d1 * k * dsin(z1), dcos(z1)],
    ]
    m2 = [
        [dcos(z2), dsin(z2) / (d2 * k)],
        [-d2 * k * dsin(z2), dcos(z2)],
    ]
    m00 = m2[0][0] * m1[0][0] + m2[0][1] * m1[1][0]
    m01 = m2[0][0] * m1[0][1] + m2[0][1] * m1[1][1]
    kap = -k / dtan(k / 2)
    value = m00 + m01 * kap
    return value.v, value.d, kap.v, kap.d


w1 = I(25) / 66
w2 = I(4) / 33
D1 = I(25)
D2 = I(1) / 4


def transfer(k):
    z1 = k * w1 / iv.sqrt(D1)
    z2 = k * w2 / iv.sqrt(D2)
    m1 = [
        [iv.cos(z1), iv.sin(z1) / (iv.sqrt(D1) * k)],
        [-iv.sqrt(D1) * k * iv.sin(z1), iv.cos(z1)],
    ]
    m2 = [
        [iv.cos(z2), iv.sin(z2) / (iv.sqrt(D2) * k)],
        [-iv.sqrt(D2) * k * iv.sin(z2), iv.cos(z2)],
    ]
    return [
        [m2[0][0] * m1[0][0] + m2[0][1] * m1[1][0],
         m2[0][0] * m1[0][1] + m2[0][1] * m1[1][1]],
        [m2[1][0] * m1[0][0] + m2[1][1] * m1[1][0],
         m2[1][0] * m1[0][1] + m2[1][1] * m1[1][1]],
    ]


def kappa_from_right_root(k):
    return -k / iv.tan(k / 2)


def crossing_characteristic(k):
    m = transfer(k)
    kap = kappa_from_right_root(k)
    return m[0][0] + m[0][1] * kap


def left_characteristic(k, kap):
    m = transfer(k)
    return m[0][0] + m[0][1] * kap


def right_characteristic(k, kap):
    return iv.cos(k / 2) + (kap / k) * iv.sin(k / 2)


def bracket_check(name, fun, bounds, kap):
    lo, hi = bounds
    flo = fun(I(lo), kap)
    fhi = fun(I(hi), kap)
    assert flo.b < 0 < fhi.a or flo.a > 0 > fhi.b, (name, flo, fhi)
    return flo, fhi


klo, khi = DATA["crossing_wave_number_bracket"]
flo = crossing_characteristic(I(klo))
fhi = crossing_characteristic(I(khi))
assert flo.b < 0 and fhi.a > 0
k_interval = I([klo, khi])
kap_interval = kappa_from_right_root(k_interval)
_, cross_derivative, _, kappa_derivative = crossing_derivative([klo, khi])
assert cross_derivative.a > 0 and kappa_derivative.a > 0

checks = {}
for name in ("L2", "L3", "L4"):
    checks[name] = bracket_check(name, left_characteristic,
                                 DATA["root_brackets"][name], kap_interval)
for name in ("R2", "R3", "R4"):
    checks[name] = bracket_check(name, right_characteristic,
                                 DATA["root_brackets"][name], kap_interval)

root_sq = {}
for name, (lo, hi) in DATA["root_brackets"].items():
    k = I([lo, hi])
    root_sq[name] = k * k

assert root_sq["R2"].b < root_sq["L2"].a
assert root_sq["L2"].b < root_sq["R3"].a
assert root_sq["R3"].b < root_sq["L3"].a
assert root_sq["R4"].b < root_sq["L4"].a

ell_l = Fraction(25, 66) / 5 + Fraction(4, 33) / Fraction(1, 2)
ell_r = Fraction(1, 2)
delta = ell_r - ell_l
assert ell_l == Fraction(7, 22)
assert delta == Fraction(2, 11)
assert 7 * delta > 1

print("INTERVAL CERTIFICATE: first left and right roots coincide.")
print("kappa crossing bracket:", kap_interval)
print("common eigenvalue bracket:", k_interval * k_interval)
print("crossing characteristic at lower endpoint:", flo)
print("crossing characteristic at upper endpoint:", fhi)
print("crossing characteristic derivative interval:", cross_derivative)
print("dkappa/dk interval:", kappa_derivative)
for name in ("R2", "L2", "R3", "L3", "R4", "L4"):
    print(name, "k^2 interval:", root_sq[name])
print("strict order: R2 < L2 < R3 < L3 and R4 < L4")
print("optical lengths:", ell_l, ell_r, "gap:", delta)
print("high-mode phase condition: 7*(ell_R-ell_L) > 1")
print("All interval assertions passed.")
print("Python:", platform.python_version())
print("mpmath:", __import__("mpmath").__version__)
print("Platform:", platform.platform())
try:
    print("CPU:", subprocess.check_output(
        ["sh", "-c", "grep -m1 'model name' /proc/cpuinfo | cut -d: -f2"],
        text=True).strip())
except Exception:
    pass
