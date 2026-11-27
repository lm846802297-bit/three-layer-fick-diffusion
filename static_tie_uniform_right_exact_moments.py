"""Exact interval certificate for both suffix moments on the static-tie benchmark.

Run:
    python static_tie_uniform_right_exact_moments.py

Only mpmath is required. The low-mode root brackets are inherited from
static_tie_uniform_right_check.py; this script rechecks their residual signs
and the phase-order conditions before evaluating the exact characteristic
functions. No high-mode truncation is used in M_gamma or M_m.
"""
from fractions import Fraction as Q
from mpmath import iv

iv.dps = 50
I = iv.mpf
pi = iv.pi

# Physical benchmark: left fast layer, left slow layer, uniform right half.
w1q, w2q = Q(25, 66), Q(4, 33)
D1q, D2q = Q(25), Q(1, 4)
kappaq = Q(4)
assert w1q + w2q == Q(1, 2)
assert w1q / D1q + w2q / D2q == Q(1, 2)
ellq = w1q / Q(5) + w2q / Q(1, 2)
assert ellq == Q(7, 22) and ellq < Q(1, 2)

w1, w2 = I(25) / 66, I(4) / 33
D1, D2, kap = I(25), I(1) / 4, I(4)
ell = w1 / iv.sqrt(D1) + w2 / iv.sqrt(D2)
assert 7 * (I(1) / 2 - ell) > 1

def left_phase(k):
    theta = k * (I(1) / 2 - w2) / 5 - iv.atan2(kap / 5, k)
    if theta > -pi / 2 and theta < pi / 2:
        j = 0
    elif theta > pi / 2 and theta < 3 * pi / 2:
        j = 1
    else:
        raise ValueError("left phase branch not certified")
    x = theta - j * pi
    return j * pi + iv.atan2(10 * iv.sin(x), iv.cos(x)) + 2 * k * w2

def right_eq(k, n):
    return k - (2 * n - 1) * pi - 2 * iv.atan2(kap, k)

def bracket(lo, hi, residual):
    assert residual(I(lo)) < 0 and residual(I(hi)) > 0
    return I([lo, hi])

# The same certified low-mode brackets used by the original right-only check.
right_k = [
    bracket("4.57785945620", "4.57785945622",
            lambda k: right_eq(k, 1)),
    bracket("10.17397018820", "10.17397018821",
            lambda k: right_eq(k, 2)),
    bracket("16.19232720644", "16.19232720646",
            lambda k: right_eq(k, 3)),
]
for n, k in enumerate(right_k, 1):
    margin = left_phase(k) - (I(n) - I(1) / 2) * pi
    assert (0 < margin < pi) if n == 1 else (-pi < margin < 0)

left_k1 = bracket(
    "3.77267640109", "3.77267640111",
    lambda k: left_phase(k) - pi / 2
)
left_k2 = bracket(
    "13.25724841026", "13.25724841029",
    lambda k: left_phase(k) - 3 * pi / 2
)

a = left_k1**2
L2 = left_k2**2
b, m, R3 = (k**2 for k in right_k)
q = 1 - m / L2
gamma = I(183)
assert a < m < gamma < R3
assert q > 0

def left_characteristic(k):
    z1 = k * w1 / iv.sqrt(D1)
    M1 = [
        [iv.cos(z1), iv.sin(z1) / (iv.sqrt(D1) * k)],
        [-iv.sqrt(D1) * k * iv.sin(z1), iv.cos(z1)],
    ]
    z2 = k * w2 / iv.sqrt(D2)
    M2 = [
        [iv.cos(z2), iv.sin(z2) / (iv.sqrt(D2) * k)],
        [-iv.sqrt(D2) * k * iv.sin(z2), iv.cos(z2)],
    ]
    v1 = [
        M1[0][0] + M1[0][1] * kap,
        M1[1][0] + M1[1][1] * kap,
    ]
    v2 = [
        M2[0][0] * v1[0] + M2[0][1] * v1[1],
        M2[1][0] * v1[0] + M2[1][1] * v1[1],
    ]
    return v2[0]

def P_left_from_k(k):
    # Static tie gives Y_L(0)=1+kappa/2=3.
    return left_characteristic(k) / 3

def P_right_from_gamma(g):
    k = iv.sqrt(g)
    return (iv.cos(k / 2) + (kap / k) * iv.sin(k / 2)) / 3

# Exact M_gamma from the two normalized Hadamard products.
P_L_gamma = P_left_from_k(iv.sqrt(gamma))
P_R_gamma = P_right_from_gamma(gamma)
M_gamma = (
    (P_L_gamma / P_R_gamma)
    * ((1 - gamma / b) * (1 - gamma / m))
    / ((1 - gamma / a) * (1 - gamma / L2))
)

# Exact M_m from the simple-pole limit at m=R_2.
P_L_m = P_left_from_k(right_k[1])
dP_R_dk = (
    -iv.sin(right_k[1] / 2) / 2
    + kap * (
        iv.cos(right_k[1] / 2) / (2 * right_k[1])
        - iv.sin(right_k[1] / 2) / (right_k[1]**2)
    )
) / 3
P_R_prime_at_minus_m = -dP_R_dk / (2 * right_k[1])
M_m = (
    P_L_m * (1 - m / b)
    / (
        m * P_R_prime_at_minus_m
        * (1 - m / a) * (1 - m / L2)
    )
)

# Insert exact moments into the lower bound for f-aS.
cutoff = iv.ln(
    M_gamma * (1 + a / (q * (m - a))) / M_m
) / (gamma - m)
assert cutoff < I("0.0171")

print("Static tie: exact resistance match R_L=R_R=1/2.")
print("Left optical length:", ellq)
for name, value in [
    ("a=L1", a),
    ("b=R1", b),
    ("m=R2", m),
    ("L2", L2),
    ("R3", R3),
    ("q=1-m/L2", q),
    ("P_L(-183)", P_L_gamma),
    ("P_R(-183)", P_R_gamma),
    ("M_183=E exp(183 T)", M_gamma),
    ("M_m=E exp(m T)", M_m),
    ("cutoff", cutoff),
]:
    print(name, value)
print("CERTIFIED: h(t)>a for every t >= 0.0171.")
