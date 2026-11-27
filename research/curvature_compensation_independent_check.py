"""Independent high-precision root solve for the Robin crossing benchmark."""
import mpmath as mp
mp.mp.dps = 70
w1 = mp.mpf(25) / 66
w2 = mp.mpf(4) / 33


def left_matrix(k):
    z1 = k * w1 / 5
    z2 = 2 * k * w2
    m1 = mp.matrix([[mp.cos(z1), mp.sin(z1) / (5*k)],
                    [-5*k*mp.sin(z1), mp.cos(z1)]])
    m2 = mp.matrix([[mp.cos(z2), 2*mp.sin(z2)/k],
                    [-k*mp.sin(z2)/2, mp.cos(z2)]])
    return m2*m1


def kappa_from_right(k):
    return -k/mp.tan(k/2)


def crossing(k):
    m = left_matrix(k)
    kap = kappa_from_right(k)
    return m[0,0] + kap*m[0,1]


k0 = mp.findroot(crossing, (mp.mpf('5.3'), mp.mpf('5.4')))
kap0 = kappa_from_right(k0)


def left(k):
    m = left_matrix(k)
    return m[0,0] + kap0*m[0,1]


def right(k):
    return mp.cos(k/2) + (kap0/k)*mp.sin(k/2)


def roots(fun, count, xmax=100):
    vals = []
    ngrid = 200000
    x0 = mp.mpf('1e-8')
    y0 = fun(x0)
    for j in range(1, ngrid + 1):
        x = mp.mpf(xmax)*j/ngrid
        y = fun(x)
        if y*y0 < 0:
            root = mp.findroot(fun, (x0, x))
            if not vals or root - vals[-1] > mp.mpf('1e-6'):
                vals.append(root)
            if len(vals) == count:
                break
        x0 = x
        y0 = y
    return vals


L = roots(left, 4)
R = roots(right, 4)
print('kappa0 =', mp.nstr(kap0, 60))
print('k0 =', mp.nstr(k0, 60))
for i, (kl, kr) in enumerate(zip(L, R), 1):
    print(f'L{i} =', mp.nstr(kl**2, 45))
    print(f'R{i} =', mp.nstr(kr**2, 45))
print('crossing derivative =', mp.nstr(mp.diff(crossing, k0), 40))
print('dkappa/dk =', mp.nstr(mp.diff(kappa_from_right, k0), 40))
