from mpmath import iv
iv.dps=40
I=iv.mpf; pi=iv.pi
z=I(".12"); v=I(".02"); u=I(2); kappa=I("5.0902282130")
ell1=(I(".5")-z)*v; ell2=z*u
def phase(k):
    d=k*ell1-iv.atan2(kappa*v,k)
    assert -pi/2<d<pi/2
    return iv.atan2((u/v)*iv.sin(d),iv.cos(d))+k*ell2
def br(lo,hi,fun,inc=True):
    sg=1 if inc else -1
    assert sg*fun(I(lo))<0<sg*fun(I(hi))
    return I([lo,hi])
def right(lo,hi,n):
    return br(lo,hi,lambda k:k-(2*n-1)*pi-2*iv.atan2(kappa,k))
a=br("4.11841382541","4.11841382544",lambda k:phase(k)-pi/2)**2
L2=br("13.52446920049","13.52446920052",lambda k:phase(k)-3*pi/2)**2
b=right("4.77605392660","4.77605392663",1)**2
m=right("10.33971330418","10.33971330421",2)**2
assert a<b and L2>m
# n>=3 follows from the standard optical-phase gap estimate.
assert 5*(I(".5")-ell1-ell2)>1
RR=1/kappa+I(".5"); RL=1/kappa+(I(".5")-z)*v*v+z*u*u
gamma=RR/RL
U=(kappa+6)/(24*(kappa+2))-1/b
r=a*U
assert gamma>1 and gamma*a/b < 1-b*U*(1-a/b)
Fplus=(1+r+iv.sqrt((1+r)**2-4*r*gamma))/2
q=1-m/L2
TL=iv.ln(q/(1-Fplus))/m
assert TL>I(".0298")
print("early-window lower bound",TL)

def W(k):
    x,y,w=k*ell1,k*ell2,k/2
    s0,s1,s2=iv.sin(x),iv.sin(y),iv.sin(w)
    a0,a1,a2=iv.cos(x),iv.cos(y),iv.cos(w)
    Y=a0+kappa*v/k*s0; Q=kappa*a0-k/v*s0
    Yp=-ell1*s0+kappa*v*(ell1*a0/k-s0/k**2)
    Qp=-kappa*ell1*s0-(s0+k*ell1*a0)/v
    L=Y*a1+u/k*Q*s1
    A=Q*a1-k/u*Y*s1
    Lp=Yp*a1-Y*ell2*s1+u/k*(Qp*s1+Q*ell2*a1)-u*Q*s1/k**2
    Ap=Qp*a1-Q*ell2*s1-(Y*s1+k*Yp*s1+k*Y*ell2*a1)/u
    R=a2+kappa/k*s2
    S=kappa*a2-k*s2
    Rp=-s2/2+kappa*(a2/(2*k)-s2/k**2)
    Sp=-kappa*s2/2-s2-k*a2/2
    return L*S+R*A,L,R,Lp*S+L*Sp+Rp*A+R*Ap
k1=br("2.64390422161","2.64390422165",lambda k:W(k)[0],False)
k2=br("4.26468719976","4.26468719979",lambda k:W(k)[0])
lam,c=k1**2,k2**2
assert lam<a and a<c<b
_,L,R,der=W(k2)
YL0=kappa*RL
YR0=kappa*RR
W0=kappa*(YL0+YR0)
Pprime=-der/(2*k2*W0)
M=(R/YR0)*(1-c/lam)/((1-c/b)*c*Pprime)
nu=gamma*(L/R)*(1-c/b)/(1-c/a)
Binf=1/b-1/c+(1/c-1/a)*nu
late=iv.ln(M*(2-c/b-c*Binf/(gamma-1)))/c
print("late-window cutoff",late)
assert late<I(".23")
# Abstract single-valley shape counterexample, signs enclosed.
assert iv.exp(I("-.2"))*(1-iv.exp(I(-1)))-2*(1-iv.exp(I("-.2")))>I(".02")
assert iv.exp(I("-4.8"))-I(".01")*(1-iv.exp(I("-4.8")))<0
assert -2*iv.exp(I(-1))-I(".01")*(1-iv.exp(I(-1)))>-1
print("generic exponential-stage four-crossing signs verified")
