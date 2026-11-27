"""Two-stage half-spectrum gate and Abel full-spectrum bridge."""
from functools import lru_cache
from fractions import Fraction as Q
import math as m
from scipy.optimize import brentq
from mpmath import iv

iv.dps=40
I=iv.mpf; pi=iv.pi
N=24
L=[".38",".12",".5"]; V=[".02","2","1"]
e1=I(L[0])*I(V[0]); e2=I(L[1])*I(V[1])
assert 5*(I(".5")-e1-e2)>1
bulk=I(L[0])*I(V[0])**2+I(L[1])*I(V[1])**2
assert bulk<I(".5")  # d positive at every kappa

def rat(x):return I(x.numerator)/x.denominator
def join(x,y):return I([x.a,y.b])
def phase(k,K,numpy=False):
    atan=m.atan2 if numpy else iv.atan2
    cos=m.cos if numpy else iv.cos
    sin=m.sin if numpy else iv.sin
    E1=.0076 if numpy else e1; E2=.24 if numpy else e2
    x=k*E1-atan(K*.02,k) if numpy else k*E1-atan(K*I(".02"),k)
    h=m.pi if numpy else pi
    assert -h/2<x<h/2
    return atan(100*sin(x),cos(x))+k*E2

def calc(k,K,cert=False,numpy=False,count_nodes=False):
    cast=float if numpy else I
    sin=m.sin if numpy else iv.sin
    cos=m.cos if numpy else iv.cos
    y=cast(1);p=K/k;norm=cast(0);nodes=0
    for j,(W,r0) in enumerate(zip(L,V)):
        w=cast(W);r=cast(r0)
        a,b,x=y,r*p,r*w*k
        C,S=cos(x),sin(x)
        if cert:
            norm+=w*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*k)
        if count_nodes:
            ph=iv.atan2(b,a)
            offs=[(x-ph-pi/2)/pi,(-ph-pi/2)/pi]
            f=[m.floor(float(v.a)) for v in offs]
            assert all(lo < v < lo+1 for lo,v in zip(f,offs))
            nodes += f[0]-f[1]
        y,p=a*C+b*S,(b*C-a*S)/r
        if j==1: mid=y
    if cert:
        assert norm>0
        return K*mid*(1-y)/(norm*k*k),nodes
    return p+K*y/k

def cert_bracket(x, fun):
    scale=10**10
    n=m.floor(x*scale)
    l,r=rat(Q(n-1,scale)),rat(Q(n+2,scale))
    fl,fr=fun(l),fun(r)
    assert (fl<0<fr) or (fr<0<fl)
    return join(l,r)

@lru_cache(None)
def roots(kq):
    K=rat(kq);kf=float(kq)
    rr=[]
    for n in [1,2,3]:
        x=brentq(lambda k:k-(2*n-1)*m.pi-2*m.atan2(kf,k),(2*n-1)*m.pi,2*n*m.pi)
        rr.append(cert_bracket(x,lambda k:k-(2*n-1)*pi-2*iv.atan2(K,k)))
    aa=[]
    for n in [1,2,3]:
        x=brentq(lambda k:phase(k,kf,True)-(n-.5)*m.pi,.001,45)
        aa.append(cert_bracket(x,lambda k:phase(k,K)-(I(n)-I(".5"))*pi))
    a=.001;f=calc(a,kf,numpy=True);found=[]
    while len(found)<N and a<300:
        b=a+.15;g=calc(b,kf,numpy=True)
        if f*g<0:
            x=brentq(lambda k:calc(k,kf,numpy=True),a,b)
            kk=cert_bracket(x,lambda k:calc(k,K))
            _,id=calc(kk,K,True,count_nodes=True)
            assert id==len(found)
            found.append(kk)
        a,f=b,g
    assert len(found)==N
    return aa,rr,found

def verify(ka,kb):
    A,R,G=roots(ka); B,S,H=roots(kb)
    K=join(rat(ka),rat(kb))
    a,l,l3=[join(x,y)**2 for x,y in zip(A,B)]; b,m2,m3=[join(x,y)**2 for x,y in zip(R,S)]
    assert a<b and l>m2 and l3>m3 and m3>m2
    RR=I(".5")+1/K; RL=bulk+1/K; gamma=RR/RL
    U=(K+6)/(24*(K+2))-1/b; r=a*U
    assert U>0 and gamma*a/b<1-b*U*(1-a/b)
    Fplus=(1+r+iv.sqrt((1+r)**2-4*r*gamma))/2
    assert Fplus<1
    q=1-m2/l
    q3=1-m3/l3
    def survival(t):
        e2=iv.exp(-m2*t);e3=iv.exp(-m3*t)
        return q*(1-q3)*e2+(1-q)*q3*e3+q*q3*(m3*e2-m2*e3)/(m3-m2)
    a0,b0=0,10000
    while b0-a0>1:
        k=(a0+b0)//2
        if survival(rat(Q(k,100000)))>1-Fplus:a0=k
        else:b0=k
    T=rat(Q(a0,100000))
    assert I(".02")<T and survival(T)>1-Fplus
    eigen=[];cs=[]
    for i,(x,y) in enumerate(zip(G,H)):
        k=join(x,y)
        if i: assert k*k>eigen[-1]
        coef,_=calc(k,K,True)
        cs.append(coef);eigen.append(k*k)
    L1,L2,L3=eigen[:3]
    assert a<L2<b and cs[1]<0 and cs[2]>0
    margin=cs[2]*L3*(L3-L1)*(L3-L2)
    margin_floor=margin.a
    for lam,cn in zip(eigen[3:],cs[3:]):
        margin+=cn*lam*(lam-L1)*(lam-L2)*iv.exp(-(lam-L3)*T)
        margin_floor=min(margin_floor,margin.a)
    # D_min*pi^2>2; Sobolev gives |lambda_n*c_n|<=2*kappa*(1+4 sqrt(lambda_n)).
    z=2*N*N
    assert z*T>I("2.5")
    j=I(N+1)/N
    ratio=j**5*iv.exp(-2*(2*N+1)*T); assert ratio<1
    tail=10*K*z**I("2.5")*iv.exp(-(z-L3)*T)/(1-ratio)
    margin_floor=min(margin_floor,(margin-tail).a)
    assert margin_floor>0
    return T, margin_floor

if __name__=="__main__":
    import sys
    LO,HI,STEP=map(int,sys.argv[1:] or [5350,6300,10])
    vals=[]; missing=[]
    todo=[(Q(k,1000), Q(min(k+STEP,HI),1000)) for k in range(LO,HI,STEP)]
    while todo:
        a,b=todo.pop()
        try:
            T,r=verify(a,b)
            vals.append((T,r))
        except AssertionError:
            if b-a <= Q(1,10000):
                missing.append((a,b));continue
            md=(a+b)/2
            todo.extend(((a,md),(md,b)))
    print("KAPPA",Q(LO,1000),Q(HI,1000),"boxes",len(vals),"missing",len(missing))
    print("all unresolved", missing[:10])
    print("glue time intervals",min(t.a for t,v in vals),max(t.b for t,v in vals))
    print("minimum certified lower bound",min(v.a for t,v in vals))
    assert not missing
