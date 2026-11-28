# Each occurrence of either independent wall rate is enclosed by the same interval.
# Equal-rate endpoint corners bracket eigenvalues by min-max monotonicity.
from mpmath import iv
iv.dps=55
I=iv.mpf
pi=iv.pi
l,n=I(5)/66,I(8)/33
from sys import argv
ends=argv[1:] or ['3.9','4.1']
assert len(ends)==2
kap=I(ends)
ka=kap
assert 0<kap<I('4.5')
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
a=root(0,'3.60','3.94')**2
b=root(1,'4.45','4.72')**2
la=root(2,'2.23','2.62')**2
kc=root(2,'3.77','4.20')
km=root(2,'7.20','7.40'); c=kc**2; lam=km**2
B=I(168)
assert la<a<c<b<24 and b<lam<B
k=iv.sqrt(B)
def phase(k):
    z=k*l-iv.atan2((ka/5),k)
    assert -pi/2<z<pi/2
    return iv.atan2(10*iv.sin(z),iv.cos(z))+k*n
assert pi/2<phase(iv.sqrt(b))<phase(k)<3*pi/2
assert iv.sqrt(lam)<3*pi
assert pi<b**I('.5')<2*pi
assert k>3*pi+2*iv.atan2(kap,k) and k<5*pi
kn=root(2,'11.40','11.75')
lam4=kn**2
assert lam<lam4<B
assert real(k)[2]<0 and real(I('7.40'))[2]>0
# Above: there are exactly two split eigenvalues below lambda_3 and three
# below B. The fourth root bracket accounts for the only possible further root.
def amp(k,j,r):
    L,R,z,d,Rp=real(k); Y=[L,R][j]; v=k**2
    return Y*2/d*(la-v)/(r-v)
K=[amp(kc,j,r) for j,r in [(1,b),(0,a)]]
lead=[-amp(km,j,r)/q for (j,r),q in zip([(1,b),(0,a)],K)]
assert K[0]>0 and K[1]>0 and lead[0].a>0
second=[amp(kn,j,r)/q for (j,r),q in zip([(1,b),(0,a)],K)]
Dsecond=(lam4-c)**2*(second[0]-second[1])
Tsecond=-(lam4-c)**2*second[1]
print('second signed limits',Dsecond,Tsecond,flush=True)
print("eta, survival amplitudes, CDF normalization",lam-c,lead,K,flush=True)

# Bound the second derivative of the normalized CDF contrast before inversion.
X=8; y0=iv.sqrt(B+X*X); N=2048
central=[I(0)]*3
contrast_cell_bound=I(0)
for i in range(N):
    x=I([i,i+1])*X/N; y=iv.sqrt(B+x*x)
    q=iv.mpc(x,y); p=iv.mpc(-B,2*x*y)
    L,R,D=data(q)
    assert abs(D).a>0
    w=iv.sqrt(1+(x/y)**2)*(2*X/N)/pi*abs(p+c)**2
    terms=[w*abs(Y/D*(p+la)/(p+r))/K[j]
           for j,(Y,r) in enumerate(((R,b),(L,a)))]
    num=R*(p+a)*K[1]-L*(p+b)*K[0]
    both=w*abs(p+la)/abs(D)*abs(num)/(abs(p+a)*abs(p+b)*K[0]*K[1])
    for j,bnd in enumerate((*terms,both)): central[j]+=bnd
    contrast_cell_bound+=min(both.b,(terms[0]+terms[1]).b)
bd=(1+iv.exp(-16*l))*(1+iv.exp(-16*n))*(1+iv.exp(-8))
den=I(33)/16*(2-bd); assert den>I('1.35')
nums=[(1+kap/I(X))*(1+iv.exp(-X))/2,
      (11/I(4))*(1+(ka/5)/X)*(1+iv.exp(-2*l*X))*(1+iv.exp(-2*n*X))]
d=[l+n,I(1)/2]
from math import comb,factorial
eta=lam-c
assert phase(iv.sqrt(170))<3*pi/2
r2=root(1,'10.04','10.30')**2
assert lam<r2<110
g2=eta*(170-r2)/((170-c)*(r2-lam))
assert (2*pi/(l+n))**2>380
assert iv.sqrt(270)>5*pi+2*iv.atan2(kap,iv.sqrt(270))
g3=eta*(380-270)/((380-c)*(270-lam))
gap=g2+g3+g2*g3
print("two-factor relative gap",gap,flush=True)
assert gap>0
base=eta**2*lead[0]*gap
Tlead=eta**2*lead[0]*(1+gap)
print("central second derivative integrals X, T, X-T:",central,flush=True)
print("joint central upper",contrast_cell_bound.b,flush=True)
for cutoff in ('.055','.06','.062','.063','.064','.065'):
    tau=I(cutoff)
    ray=[2/pi/ K[j]*nums[j]/den*iv.exp(-d[j]*X)*(1+I(24)/X**2)**3*
         sum(comb(4,k)*(X+y0)**(4-k)*factorial(k)/(2*y0*tau)**(k+1) for k in range(5)) for j in (0,1)]
    Econtrast=contrast_cell_bound+sum(ray)
    qratio=iv.exp(-(B-lam)*tau)
    e2=iv.exp(-(lam4-lam)*tau)
    mcontrast=base-Econtrast*qratio+min(I(0),Dsecond.a)*e2
    mt=Tlead-(central[1]+ray[1])*qratio+min(I(0),Tsecond.a)*e2
    print(cutoff,"scaled lower bounds X-T and -T:",mcontrast.a,mt.a,flush=True)

    if cutoff=='.062': assert mcontrast.a>66 and mt.a>1261
