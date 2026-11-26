"""Analytic envelope for first-only inversion; interval check at one point."""
from first_inversion_residual_band import f,g,bracket,M,I

K_s,z_s='6.2','.12'
K,z=I.mpf(K_s),I.mpf(z_s)
v=I.mpf('.02'); h=I.mpf('.5')-z
def phase(q):
    t=q*h*v-I.atan2(K*v,q)
    c,s=I.cos(t),I.sin(t)
    return t+I.atan2(99*c*s,c*c+100*s*s)+2*z*q
q=bracket(K_s,z_s,'full',2,3)
ql=bracket(K_s,z_s,'left',13,14)
qb=bracket(K_s,z_s,'right',4,5)
qm=bracket(K_s,z_s,'right',9.5,11)
l,b,m,lam=ql**2,qb**2,qm**2,q**2
assert q<I.pi and g(q,K,z,I,'left')>0 and g(q,K,z,I,'right')>0
assert all(left>0 and right>0 and x<I.pi for left,right,x in f(q,K,z,I,'check'))
assert qb>I.pi and qb<2*I.pi and qm>3*I.pi and qm<4*I.pi
assert phase(qb)>I.pi/2 and phase(qm)<3*I.pi/2
assert I.pi<phase(ql)<2*I.pi and l>m
assert 5*(I.mpf('.5')-h*v-2*z)>1
rho=(h*v*v+4*z+1/K)/(I.mpf('.5')+1/K)
R=g(q,K,z,I,'left')/g(q,K,z,I,'right')/rho
old=(rho*b-lam)/(rho*(b-lam))
assert rho<1 and R<old
A=b*(1-rho)/lam
B=(b-lam)*(rho*R-1)/lam
C=(b-(b-lam)*R*l/(l-lam))/lam
E=rho*C*(1-m/l)
gamma=m/lam
x=(-B/(E*gamma))**(1/(gamma-1))
assert B<0 and C>0 and E>0 and gamma>1 and 0<x<1
val=A+B*x+E*x**gamma
assert val>I.mpf('.0073')
for name,z in [('rho',rho),('spectral ratio',R),('old lower bound',old),('minimizer x',x),('minimum',val)]:
    print(name,M.nstr(M.mpf(z.mid),15),'+/- at most',M.nstr(M.mpf(z.delta),3))
print('PASS no negative valley for any t > 0, despite failing prior test')
