"""Interval certificate that the Robin delay density is not bell-shaped.

Benchmark:
  Left half: D=25 on [0,25/66], D=1/4 on [25/66,1/2].
  Right half: D=1 on [1/2,1].
  Equal Robin coefficient kappa=4 at the outer wall and midpoint Dirichlet.

The certificate proves
  sqrt(L_21) > sqrt(R_33),
  sqrt(L_22) < sqrt(R_34),
which makes the spectral interval multiplicity m(u) satisfy
  m(u_+) >= 13 and m(u_-) <= 11 for u_+ < u_-.
Thus m(u)-12 changes sign at least twice, violating the
Kwaśnicki--Simon level-crossing condition.

Only mpmath interval arithmetic is used. The phase residuals locate and index
roots; independent transfer-matrix characteristic residuals check the same
endpoint sign changes.
"""
from fractions import Fraction as Q
import platform
import sys
import mpmath
from mpmath import iv

iv.dps = 60
I = iv.mpf
pi = iv.pi

# Exact physical input, retained as rational data.
w1q, w2q = Q(25, 66), Q(4, 33)
D1q, D2q, DRq = Q(25), Q(1, 4), Q(1)
kappaq = Q(4)
ellLq = w1q / Q(5) + w2q / Q(1, 2)
ellRq = Q(1, 2)
assert w1q + w2q == Q(1, 2)
assert ellLq == Q(7, 22)
assert ellRq == Q(1, 2)
assert ellLq < ellRq

w1, w2 = I(25) / 66, I(4) / 33
kap = I(4)
ell1, ell2 = I(5) / 66, I(8) / 33
impedance_ratio = I(10)  # sqrt(D1/D2)

# Continuous left Prüfer/phase representation. On the brackets below,
# theta lies in (5*pi-pi/2, 5*pi+pi/2), so the branch j=5 is fixed.
def left_theta(k):
    return k * ell1 - iv.atan2(I(4) / 5, k)

def left_phase_residual(k, n, j=5):
    theta = left_theta(k)
    x = theta - j * pi
    return j * pi + iv.atan2(impedance_ratio * iv.sin(x), iv.cos(x)) + k * ell2 - (I(n) - I(1) / 2) * pi

# Right phase equation: k/2 - arctan(kappa/k) = (n-1/2)*pi.
def right_phase_residual(k, n):
    return k / 2 - iv.atan2(kap, k) - (I(n) - I(1) / 2) * pi

# Independent direct transfer-matrix characteristic functions.
def left_characteristic(k):
    z1 = k * w1 / 5
    z2 = 2 * k * w2
    y1 = iv.cos(z1) + kap * iv.sin(z1) / (5 * k)
    p1 = -5 * k * iv.sin(z1) + kap * iv.cos(z1)
    return iv.cos(z2) * y1 + 2 * iv.sin(z2) * p1 / k

def right_characteristic(k):
    return iv.cos(k / 2) + kap * iv.sin(k / 2) / k

def interval(lo, hi):
    return I([lo, hi])

def assert_negative(x, label):
    assert x.b < 0, (label, x)

def assert_positive(x, label):
    assert x.a > 0, (label, x)

def certify_root(name, n, lo, hi, phase_fun, char_fun, branch=None):
    lo_i, hi_i = I(lo), I(hi)
    rlo = phase_fun(lo_i, n) if branch is None else phase_fun(lo_i, n, branch)
    rhi = phase_fun(hi_i, n) if branch is None else phase_fun(hi_i, n, branch)
    assert_negative(rlo, name + " phase lower endpoint")
    assert_positive(rhi, name + " phase upper endpoint")
    clo, chi = char_fun(lo_i), char_fun(hi_i)
    # Direct transfer characteristic must also switch sign.
    assert clo.b < 0 < chi.a or chi.b < 0 < clo.a, (name, clo, chi)
    return interval(lo, hi)

# Root brackets are narrow, decimal-rational intervals.
kL21 = certify_root(
    "L21", 21, "205.1617247434", "205.1617247438",
    left_phase_residual, left_characteristic, branch=5
)
kL22 = certify_root(
    "L22", 22, "209.5748600934", "209.5748600939",
    left_phase_residual, left_characteristic, branch=5
)
kR33 = certify_root(
    "R33", 33, "204.2426865656", "204.2426865661",
    right_phase_residual, right_characteristic
)
kR34 = certify_root(
    "R34", 34, "210.5247035095", "210.5247035100",
    right_phase_residual, right_characteristic
)

# Verify the left branch is fixed throughout both root brackets.
for name, k in [("L21", kL21), ("L22", kL22)]:
    theta = left_theta(k)
    x = theta - 5 * pi
    assert x.a > -pi / 2 and x.b < pi / 2, (name, theta, x)

# The phase maps are strictly increasing:
# theta'(k)=ell1+(4/5)/(k^2+(4/5)^2)>0;
# d/dx atan2(10 sin x, cos x)=10/(cos^2 x+100 sin^2 x)>0;
# hence left_phase'(k)>0. The right phase derivative is
# 1/2+4/(k^2+16)>0. Therefore the target phase levels index the roots.

# Positive and negative shifted spectral gaps.
d21_12 = kL21 - kR33
d22_12 = kL22 - kR34
assert_positive(d21_12, "sqrt(L21)-sqrt(R33)")
assert_negative(d22_12, "sqrt(L22)-sqrt(R34)")

# The same signs in the unshifted spectral rates (subtracting c=lambda_2
# cancels from each comparison).
rate_gap_positive = kL21**2 - kR33**2  # L_21 - R_33 > 0
rate_gap_negative = kR34**2 - kL22**2  # R_34 - L_22 > 0
assert_positive(rate_gap_positive, "L21-R33")
assert_positive(rate_gap_negative, "R34-L22")

# Spectral ordering gives rho_33 < sigma_21 < sigma_22 < rho_34.
# Choose u+ in (rho_33,sigma_21), u- in (sigma_22,rho_34).
# At u+, all 13 intervals j=21,...,33 are active, so m(u+)>=13.
# At u-, j<=22 have ended and j>=34 have not started, so m(u-)<=11.
# Thus m(u+)-12>0>m(u-)-12, with u+<u-.

print("ROBIN_DELAY_BELLSHAPE_LEVEL12_OBSTRUCTION: PASS")
print("geometry: w1=25/66, D1=25; w2=4/33, D2=1/4; right D=1; kappa=4")
print("optical lengths: ell_L=", ellLq, "ell_R=", ellRq)
print("runtime: Python", sys.version.split()[0], "mpmath", mpmath.__version__, platform.platform())
for name, value in [("sqrt(L21)", kL21), ("sqrt(L22)", kL22),
                    ("sqrt(R33)", kR33), ("sqrt(R34)", kR34),
                    ("d_21,12", d21_12), ("d_22,12", d22_12),
                    ("L21-R33", rate_gap_positive), ("R34-L22", rate_gap_negative)]:
    print(name, value)
print("spectral order: rho_33 < sigma_21 < sigma_22 < rho_34")
print("multiplicity: m(u_plus)>=13; m(u_minus)<=11; u_plus<u_minus")
print("level k=12: m(u_plus)-12>0>m(u_minus)-12")
print("CONCLUSION: the Kwaśnicki-Simon level-crossing condition fails; f_Z is not bell-shaped.")
