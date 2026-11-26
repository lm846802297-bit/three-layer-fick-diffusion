"""Interval checks of spectral ordering and smoothed coefficient signs."""
import mpmath as m
import numpy as np
from scan_critical import mode_factory
M,I=m.mp,m.iv
M.dps=73;I.dps=65
P=I.mpf
ka,kb='5.0902282130','5.0902282130'
def show(x):return M.nstr(x,70)
def calc(q,K,A,full=False):
    M0=A.mpf
    p=K/q;y=M0(1);norm=M0(0);N=[]
    for j,(ww,rr) in enumerate((('.38','.02'),('.12','2'),('.5','1'))):
        w,r=M0(ww),M0(rr)
        a,b,x=y,r*p,q*w*r; C,S=A.cos(x),A.sin(x)
        if full:
            norm+=w*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*q)
            th=A.atan2(b,a);N.append(((x-th-A.pi/2)/A.pi,(-th-A.pi/2)/A.pi))
        y,p=a*C+b*S,(b*C-a*S)/r
        if j==1:mid=y
    return (K*mid*(1-y)/(norm*q*q),N) if full else p+K*y/q
def eigbr(guess,K,eps='1e-47'):
    q=M.findroot(lambda q:calc(q,K,M),(guess-.001,guess+.001))
    a,b=[show(q+d*M.mpf(eps)) for d in (-1,1)]
    assert calc(P(a),P(show(K)),I)*calc(P(b),P(show(K)),I)<0
    return a,b
def cocoe(aa,bb,N):
    a,b=map(M.mpf,(aa,bb));K=P([aa,bb]); mid=(a+b)/2
    _,_,modes=mode_factory(2,.02,.12,False)
    L=[];C=[]
    def fl(q):
        u,v=M.mpf(q.a),M.mpf(q.b);assert M.floor(u)==M.floor(v) and u>M.floor(u)
        return int(M.floor(u))
    for g0 in np.sqrt(modes(float(mid),N)[0]):
        g=M.mpf(float(g0));lo,e=eigbr(g,a);f,hi=eigbr(g,b)
        lo,hi=P(lo),P(hi);q=P([M.mpf(lo.a),M.mpf(hi.b)])
        c,ns=calc(q,K,I,True)
        assert sum(fl(x)-fl(y) for x,y in ns)==len(L)
        assert not L or q*q>L[-1]
        L.append(q*q);C.append(c)
    assert len(L)==N
    return L,C,K
L,C,K=cocoe(ka,kb,20)
def rk(q,A):
    return A.cos(q/2)+A.mpf(ka)*A.sin(q/2)/q
R=[]
from scipy.optimize import brentq
import math
for i in range(1,4):
    x=brentq(lambda q:math.cos(q/2)+float(ka)*math.sin(q/2)/q,(2*i-1)*math.pi,2*i*math.pi)
    x=M.findroot(lambda q:rk(q,M),(x-.002,x+.002))
    a,b=[P(show(x+d*M.mpf('1e-44'))) for d in (-1,1)]
    assert rk(a,I)*rk(b,I)<0 and a>(2*i-1)*I.pi and b<2*i*I.pi
    R.append(P([M.mpf(a.a),M.mpf(b.b)])**2)
assert L[1]<R[0]<L[2]<R[1]<L[3]<R[2]
print('PASS lambda2 < muR1 < lambda3 < muR2 < lambda4 < muR3')
print('full 2-4, right 1-3',*[M.nstr(M.mpf(x.mid),12) for x in L[1:4]+R])
assert L[-1]-L[1]>2700
for m,grid in [(1,[10,100,260,700,1100]),(2,[10,100,260,700,1100,1700])]:
    vals=[]
    for j,x in enumerate(grid):
        v=P(0)
        for l,c in zip(L[1:],C[1:]):
            d=x-(l-L[1])
            assert d<0 or d>0
            if d>0:v+=c*l*(L[0]-l)*d**m
        assert (-1)**j*v>0
        vals.append(M.nstr(M.mpf(v.mid),12))
    print('PASS Riesz order',m,'alternating signs at',grid,':',vals)
