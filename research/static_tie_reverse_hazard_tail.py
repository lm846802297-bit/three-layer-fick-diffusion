from mpmath import iv
iv.dps=55
I=iv.mpf
pi=iv.pi
l,n=I(5)/66,I(8)/33
def data(q):
    def pair(x):
        w=iv.exp(x); z=1/w
        return (w+z)/2,(w-z)/2
    C,S=pair(l*q); c,s=pair(n*q); U,V=pair(q/2)
    y=C+I('.8')/q*S; z=5*S+4/q*C
    L=c*y+2*s*z; BL=s*y/2+c*z
    R=U+4/q*V; BR=V+4/q*U
    return L,R,L*BR+R*BL
def real(k):
    C,S=iv.cos(l*k),iv.sin(l*k); c,s=iv.cos(n*k),iv.sin(n*k); U,V=iv.cos(k/2),iv.sin(k/2)
    y=C+I('.8')/k*S; z=5*S-4/k*C
    yp=-l*S+I('.8')*(l*C/k-S/k**2)
    zp=5*l*C+4*C/k**2+4*l*S/k
    L=c*y-2*s*z; BL=s*y/2+c*z
    Lp=-n*s*y+c*yp-2*n*c*z-2*s*zp
    BLp=n*c*y/2+s*yp/2-n*s*z+c*zp
    R=U+4/k*V; BR=V-4/k*U
    Rp=-V/2+2*U/k-4*V/k**2
    BRp=U/2+2*V/k+4*U/k**2
    return L,R,L*BR+R*BL,Lp*BR+L*BRp+Rp*BL+R*BLp
def root(j,x,y):
    lo,hi=I(x),I(y)
    sg=1 if real(lo)[j]>0 else -1
    assert sg*real(lo)[j]>0>sg*real(hi)[j]
    for _ in range(90):
        mid=(lo+hi)/2; z=sg*real(mid)[j]
        if z.a>0:lo=mid
        elif z.b<0:hi=mid
        else:raise ArithmeticError("wide root")
    return I([lo.a,hi.b])
a=root(0,'3.77267','3.77269')**2
b=root(1,'4.57','4.59')**2
la=root(2,'2.41899','2.41900')**2
kc=root(2,'3.9','4')
km=root(2,'7.29','7.30'); c=kc**2; lam=km**2
from sys import argv
B=I(argv[1]) if len(argv)>1 else I(110)
t0=I(argv[2]) if len(argv)>2 else I('.12')
assert la<a<c<b<lam<B
k=iv.sqrt(B)
def phase(k):
    z=k*l-iv.atan2(I('.8'),k)
    assert -pi/2<z<pi/2
    return iv.atan2(10*iv.sin(z),iv.cos(z))+k*n
assert pi/2<phase(iv.sqrt(b))<phase(k)<3*pi/2
assert iv.sqrt(lam)<3*pi
assert pi<b**I('.5')<2*pi
assert k>3*pi+2*iv.atan2(4,k) and k<5*pi
assert real(k)[2]>0 and real(I('7.30'))[2]>0
# Above: there are exactly two split eigenvalues below lambda_3 and three
# below B; hence at most one further full eigenvalue up to B. Signs exclude it.
def amp(k,j,r):
    L,R,z,d=real(k); Y=[L,R][j]; v=k**2
    return Y*2/d*(la-v)/(r-v)
K=[amp(kc,j,r) for j,r in [(1,b),(0,a)]]
lead=[-amp(km,j,r)/q for (j,r),q in zip([(1,b),(0,a)],K)]
assert K[0]>0 and K[1]>0 and lead[1]>lead[0]>0
print("eta, survival amplitudes, CDF normalization",lam-c,lead,K,flush=True)
# Shifted contour: q=x+i sqrt(B+x^2), 0<=x<=8, and then q=8+iy.
# Bound L1 of the central segment using interval rectangular enclosures.
X=8; y0=iv.sqrt(B+X*X); E=[[I(0)]*3 for i in range(2)]
from fractions import Fraction
N=2048
for i in range(N):
    x=I([i,i+1])*X/N; y=iv.sqrt(B+x*x)
    q=iv.mpc(x,y); p=iv.mpc(-B,2*x*y)
    L,R,D=data(q)
    assert abs(D).a>0,(i,D)
    for s,(Y,r) in enumerate(((R,b),(L,a))):
        w=abs(Y/D*(p+la)/(p+r))*iv.sqrt(1+(x/y)**2)*(2*X/N)/pi/K[s]
        for j in range(3):E[s][j]+=w*abs(p+c)**j
print("central contour masses",E,flush=True)
# Sharp tail: on x=8 absolute ratio estimates from the same 8-exponential
# expansion as the finite-time certificate. Separate algebraic degree per jet.
bd=(1+iv.exp(-16*l))*(1+iv.exp(-16*n))*(1+iv.exp(-8))
den=I(33)/16*(2-bd); assert den>I('1.35')
nums=[(1+4/I(X))*(1+iv.exp(-X))/2,
      (11/I(4))*(1+I('.8')/X)*(1+iv.exp(-2*l*X))*(1+iv.exp(-2*n*X))]
d=[l+n,I(1)/2]
from math import comb,factorial
for s in range(2):
    for j in range(3):
        g=2/pi/ K[s]*nums[s]/den*iv.exp(-d[s]*X)*(1+I(21)/X**2)**(j+1)
        # Put u=y-y0. Use exp(-2*y0*t0*u), dropping exp(-t0*u^2).
        E[s][j]+=g*sum(comb(2*j,k)*(X+y0)**(2*j-k)*factorial(k)/(2*y0*t0)**(k+1) for k in range(2*j+1))
eta=lam-c; rate=B-c
s0=iv.exp(-eta*t0); e=iv.exp(-(B-lam)*t0)
# At all t>=t0 scale f' by exp(eta*t), with denominator F in intervals.
# Using positive densities, dropping -f_T^2/F_T^2 is safe in an upper bound.
Fs=[I([max(I(0),(1-(a0+E[s][0]*e)*s0).a), I(1)]) for s,a0 in enumerate(lead)]
assert all(F.a>0 for F in Fs)
ds=[-eta**2*a0+I([-1,1])*E[s][2]*e for s,a0 in enumerate(lead)]
fXmax=eta*lead[0]+E[0][1]*e
upper=(ds[1]/Fs[1]-ds[0]/Fs[0]+s0*(fXmax/Fs[0])**2).b
print("tail-inclusive masses",E,flush=True)
print("scaled upper bound for all t >= ",t0,upper)
assert upper<0
