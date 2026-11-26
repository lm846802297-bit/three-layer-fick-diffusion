from math import pi
from scipy.optimize import brentq
import mpmath as m
m.mp.dps=65; m.iv.dps=55; M,I=m.mp,m.iv
r,u,z,K='.0529','2.8464','.12717','84.737'
def vals(A):return [A.mpf(v) for v in (r,u,z,K)]
def yR(q,A):
    _,_,_,k=vals(A)
    return A.cos(q/2)+k*A.sin(q/2)/q
def yL(q,A):
    a,b,z,k=vals(A)
    x,y=q*(A.mpf('.5')-z)*a,q*z*b
    return (A.cos(x)+a*k*A.sin(x)/q)*A.cos(y)+b*(-A.sin(x)/a+k*A.cos(x)/q)*A.sin(y)
def phase(q,A):
    a,b,z,k=vals(A); r=b/a
    th=q*(A.mpf('.5')-z)*a-A.atan2(a*k,q)
    c,s=A.cos(th),A.sin(th)
    return th+A.atan2((r-1)*c*s,c*c+r*s*s)+q*z*b
def isolate(f,guess):
    k=M.findroot(lambda q:f(q,M),guess)
    a,b=[M.nstr(k+d*M.mpf('1e-35'),60) for d in (-1,1)]
    assert f(I.mpf(a),I)*f(I.mpf(b),I)<0
    return I.mpf([a,b])
a,b,w,k=vals(I)
gap=I.mpf('.5')-a*(I.mpf('.5')-w)-b*w
assert gap>0 and 11*gap>1
rho=((I.mpf('.5')-w)*a*a+w*b*b+1/k)/(I.mpf('.5')+1/k)
assert rho>1
invs=[]; left={};right={}
for n in range(1,6):
    rr=isolate(yR,((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    assert rr>(2*n-1)*I.pi and rr<2*n*I.pi
    th=(I.mpf(n)-I.mpf('.5'))*I.pi
    delta=phase(rr,I)-th
    if delta>0:invs.append(n)
    else:assert delta<0
    if n in [2,3]:
        q=brentq(lambda q:float(phase(M.mpf(q),M)-(n-.5)*M.pi),.001,50)
        ll=isolate(yL,(q-.001,q+.001))
        assert abs(phase(ll,I)-th)<I.mpf('1e-32')
        left[n]=ll**2;right[n]=rr**2
assert invs==[3]
margin=(left[2]-right[2])-(right[3]-left[3])
assert margin>39
print('PASS only n=3 inverted, n>=6 from optical bound')
print('PASS smaller left optical length, greater final right capture')
print('muL2',left[2],'muR2',right[2])
print('muL3',left[3],'muR3',right[3])
print('positive compensation margin',margin)
