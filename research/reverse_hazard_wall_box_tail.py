# Each occurrence of either independent wall rate is enclosed by the same interval.
# Equal-rate endpoint corners bracket eigenvalues by min-max monotonicity.
from mpmath import iv
iv.dps=55
I=iv.mpf
pi=iv.pi
l,n=I(5)/66,I(8)/33
from sys import argv
ends=argv[1:] or ['3.99999','4.00001']
assert len(ends)==2
kap=I(ends)
ka=kap
assert 0<kap<I('4.01')
def data(q):
    def pair(x):
        w=iv.exp(x); z=1/w
        return (w+z)/2,(w-z)/2
    C,S=pair(l*q); c,s=pair(n*q); U,V=pair(q/2)
    y=C+(ka/5)/q*S; z=5*S+ka/q*C
    L=c*y+2*s*z; BL=s*y/2+c*z
    R=U+ka/q*V; BR=V+ka/q*U
    return L,R,L*BR+R*BL
def real(k):
    C,S=iv.cos(l*k),iv.sin(l*k); c,s=iv.cos(n*k),iv.sin(n*k); U,V=iv.cos(k/2),iv.sin(k/2)
    y=C+(ka/5)/k*S; z=5*S-ka/k*C
    yp=-l*S+(ka/5)*(l*C/k-S/k**2)
    zp=5*l*C+ka*C/k**2+ka*l*S/k
    L=c*y-2*s*z; BL=s*y/2+c*z
    Lp=-n*s*y+c*yp-2*n*c*z-2*s*zp
    BLp=n*c*y/2+s*yp/2-n*s*z+c*zp
    R=U+ka/k*V; BR=V-ka/k*U
    Rp=-V/2+ka/2*U/k-ka*V/k**2
    BRp=U/2+ka/2*V/k+ka*U/k**2
    return L,R,L*BR+R*BL,Lp*BR+L*BRp+Rp*BL+R*BLp,Rp
def root0(j,x,y):
    lo,hi=I(x),I(y)
    sg=1 if real(lo)[j]>0 else -1
    assert sg*real(lo)[j]>0>sg*real(hi)[j]
    for _ in range(90):
        mid=(lo+hi)/2; z=sg*real(mid)[j]
        if z.a>0:lo=mid
        elif z.b<0:hi=mid
        else:raise ArithmeticError("wide root")
    return I([lo.a,hi.b])
def root(j,x,y):
    global ka
    ka=kap
    assert (real(I(x))[j]*real(I(y))[j]).b<0
    ka=I(ends[0]); z0=root0(j,x,y)
    ka=I(ends[1]); z1=root0(j,x,y)
    ka=kap
    return I([z0.a,z1.b])
a=root(0,'3.77','3.78')**2
b=root(1,'4.57','4.59')**2
la=root(2,'2.41','2.43')**2
kc=root(2,'3.9','4')
km=root(2,'7.29','7.30'); c=kc**2; lam=km**2
B=I(110)
t0=I('.035')
assert la<a<c<b<21 and b<lam<B
k=iv.sqrt(B)
def phase(k):
    z=k*l-iv.atan2((ka/5),k)
    assert -pi/2<z<pi/2
    return iv.atan2(10*iv.sin(z),iv.cos(z))+k*n
assert pi/2<phase(iv.sqrt(b))<phase(k)<3*pi/2
assert iv.sqrt(lam)<3*pi
assert pi<b**I('.5')<2*pi
assert k>3*pi+2*iv.atan2(kap,k) and k<5*pi
assert real(k)[2]>0 and real(I('7.30'))[2]>0
# Above: there are exactly two split eigenvalues below lambda_3 and three
# below B; hence at most one further full eigenvalue up to B. Signs exclude it.
def amp(k,j,r):
    L,R,z,d,Rp=real(k); Y=[L,R][j]; v=k**2
    return Y*2/d*(la-v)/(r-v)
K=[amp(kc,j,r) for j,r in [(1,b),(0,a)]]
lead=[-amp(km,j,r)/q for (j,r),q in zip([(1,b),(0,a)],K)]
assert K[0]>0 and K[1]>0 and lead[1]>lead[0]>0
print("eta, survival amplitudes, CDF normalization",lam-c,lead,K,flush=True)
# Shifted contour: q=x+i sqrt(B+x^2), 0<=x<=8, and then q=8+iy.
# Bound L1 of the central segment using interval rectangular enclosures.
X=8; y0=iv.sqrt(B+X*X); E=[[I(0)]*4 for i in range(2)]
from fractions import Fraction
N=2048
for i in range(N):
    x=I([i,i+1])*X/N; y=iv.sqrt(B+x*x)
    q=iv.mpc(x,y); p=iv.mpc(-B,2*x*y)
    L,R,D=data(q)
    assert abs(D).a>0,(i,D)
    for s,(Y,r) in enumerate(((R,b),(L,a))):
        w=abs(Y/D*(p+la)/(p+r))*iv.sqrt(1+(x/y)**2)*(2*X/N)/pi/K[s]
        for j in range(4):E[s][j]+=w*abs(p+c)**j
print("central contour masses",E,flush=True)
# Sharp tail: on x=8 absolute ratio estimates from the same 8-exponential
# expansion as the finite-time certificate. Separate algebraic degree per jet.
bd=(1+iv.exp(-16*l))*(1+iv.exp(-16*n))*(1+iv.exp(-8))
den=I(33)/16*(2-bd); assert den>I('1.35')
nums=[(1+kap/I(X))*(1+iv.exp(-X))/2,
      (11/I(4))*(1+(ka/5)/X)*(1+iv.exp(-2*l*X))*(1+iv.exp(-2*n*X))]
d=[l+n,I(1)/2]
from math import comb,factorial
for s in range(2):
    for j in range(4):
        g=2/pi/ K[s]*nums[s]/den*iv.exp(-d[s]*X)*(1+I(21)/X**2)**(j+1)
        # Put u=y-y0. Use exp(-2*y0*t0*u), dropping exp(-t0*u^2).
        E[s][j]+=g*sum(comb(2*j,k)*(X+y0)**(2*j-k)*factorial(k)/(2*y0*t0)**(k+1) for k in range(2*j+1))
eta=lam-c; rate=B-c
s0=iv.exp(-eta*t0); e=iv.exp(-(B-lam)*t0)


t=I('.10')
ratios=[(E[0][j]*iv.exp(-(B-lam)*t)/(eta**j*lead[0])).b for j in (2,3)]
print('relative remainder bounds from t=.10:',*ratios)
assert ratios[0]<1 and ratios[1]<1

# Density of Z, with the cancelled pole at -R_1 removed.
r2=root(1,'10.17','10.18')
Bz=I(220)
assert 3*pi<r2<iv.sqrt(Bz)<5*pi
rho=r2**2-c
L,R,z,d,Rp=real(r2)
leadZ=-2*r2*L*(b-r2**2)/(a-r2**2)/Rp
assert leadZ>0 and 85<rho<95
fix=I('.07')
y0=iv.sqrt(Bz+X*X)
mass=[I(0),I(0)]
for i in range(N):
    x=I([i,i+1])*X/N; y=iv.sqrt(Bz+x*x)
    q=iv.mpc(x,y); p=iv.mpc(-Bz,2*x*y)
    L,R,D=data(q)
    assert abs(R).a>0
    w=abs(q*L/R*(p+b)/(p+a))*iv.sqrt(1+(x/y)**2)*(2*X/N)/pi
    mass[0]+=w*abs(p+c+85); mass[1]+=w*abs(p+c+95)
# On the external ray Q=8, use |p+c+u|<=|q|^2+116
# for u=85 or 95, since c<21. Bound it by (1+116/64)|q|^2.
qfac=1+I(116)/64
# On Re q=8 the upper half trace constants are as above; the
# right trace has lower factor (1-exp(-8))/2.
zpref=2/pi*nums[1]*2/(1-iv.exp(-8))*iv.exp(-(I(2)/11)*X)*(1+I(21)/64)*qfac
out=zpref*sum(comb(3,k)*(X+y0)**(3-k)*factorial(k)/(2*y0*fix)**(k+1) for k in range(4))
scaled=[(mm+out)*iv.exp(-(Bz-r2**2)*fix) for mm in mass]
margins=[(rho-85)*leadZ,(95-rho)*leadZ]
ratio=[(ee/mm).b for ee,mm in zip(scaled,margins)]
print('Z leading rate, residue and relative score remainders:',rho,leadZ,ratio,flush=True)
assert max(ratio)<1
