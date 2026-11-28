from mpmath import iv
from sys import argv
iv.dps=55
I=iv.mpf
def data(q,ctx):
    S=lambda x:(ctx.exp(x)+ctx.exp(-x))/2
    H=lambda x:(ctx.exp(x)-ctx.exp(-x))/2
    l=ctx.mpf(5)/66; n=ctx.mpf(8)/33
    c,s,C,t,U,V=S(l*q),H(l*q),S(n*q),H(n*q),S(q/2),H(q/2)
    y=c+ctx.mpf(4)/5/q*s; z=5*s+4/q*c
    L=C*y+2*t*z; BL=t*y/2+C*z
    R=U+4/q*V; BR=V+4/q*U
    return L,R,L*BR+R*BL
def spec(x):return data(iv.mpc(0,x),iv)[2].imag
def rs(x):return data(iv.mpc(0,x),iv)[1].real
def ls(x):return data(iv.mpc(0,x),iv)[0].real
def root(fun,x,y):
    lo,hi=I(x),I(y)
    z=fun(lo); z1=fun(hi)
    sg=1 if z.a>0 else -1
    assert (sg*z).a>0 and (sg*z1).b<0,(x,y,z,z1)
    for j in range(100):
        mid=(lo+hi)/2; v=sg*fun(mid)
        if v.a>0:lo=mid
        elif v.b<0:hi=mid
        else:raise ArithmeticError("uncertain bracket")
    return I([lo.a,hi.b])**2

a=root(ls,'3.77267640109','3.77267640111'); b=root(rs,'4.5','4.6')
la=root(spec,'2.41899445303','2.41899445305'); c=root(spec,'3.9','4')
assert la<a<c<b


# Exact constants in the short-time analytic proof (see accompanying note).
from fractions import Fraction as Fr
dx,dt=Fr(7,22),Fr(1,2)
t0=Fr(3,2000)
assert dx/(2*t0)>100 and 64*t0<dx**2
assert (dt**2-dx**2)/2-540*t0**2*(1/dx+1/dt)-6400*t0**3*(1/dx**2+1/dt**2)>Fr(67,1000)

# python script lo hi Q_X Q_T dt, all times exact decimal inputs
from fractions import Fraction
from math import comb, factorial
lo,hi,qx,qt,step=argv[1:]
T0,T1,st=map(I,(lo,hi,step)); N=(Fraction(hi)-Fraction(lo))/Fraction(step)
assert N.denominator==1
N=int(N)
w=st/2; DEG=6; late=min(int(qx),int(qt))<30
A=8 if late else 10; h=I(1)/(4 if late else 2)
XMIN=min(int(qx),int(qt))-A if late else 20
assert XMIN>=8
l=I(5)/66; n=I(8)/33
# The 8-term exponential expansion of D gives this lower-bound criterion
# throughout Re(q)>=20. The remaining two bounds verify the numerator,
# rational factor, and maximum derivative-weight constants.
B=(1+iv.exp(-2*XMIN*l))*(1+iv.exp(-2*XMIN*n))*(1+iv.exp(-XMIN))
Db=I('1.35') if late else I('1.9')
assert B<I('1.34')
assert I(33)/16*(2-B)>Db
assert (11+I('8.8')/XMIN)/Db<9
assert 21>max(a.b,b.b,c.b,la.b)
assert 9/iv.pi*(1+I(21)/XMIN**2)**3<10
mom=[iv.sqrt(iv.pi/T0)/2,1/(2*T0)]
for j in range(2,5):mom.append((j-1)*mom[j-2]/(2*T0))
def setup(qx,Yr):
    Q=int(qx); assert Q-A>=XMIN
    delta=l+n if Yr else I(1)/2
    mx=Q+A
    U=20*iv.exp(T1*mx**2-delta*(Q-A))*sum(comb(4,j)*mx**(4-j)*mom[j] for j in range(5))
    rr=iv.exp(-2*iv.pi*A/h)
    alias=2*U*rr/(1-rr)
    M=400 if late else 500; y=(M+1)*h
    rr=((Q+y+h)/(Q+y))**4*iv.exp(-T0*(2*y*h+h*h))
    assert rr<1
    tail=20*h*iv.exp(T1*Q**2-delta*Q)*(Q+y)**4*iv.exp(-T0*y*y)/(1-rr)
    E=(alias+tail).b
    err=I([-E,E])
    data_rows=[]
    for j in range(M+1):
        q=iv.mpc(Q,j*h); p=q*q
        L,R,D=data(q,iv)
        Y,r=(R,b) if Yr else (L,a)
        z=(1 if not j else 2)*h/iv.pi*Y/D*(p+la)/(p+r)
        z*=iv.exp(p*(T0+w))
        rho=iv.exp(p*st)
        # Remainder at the midpoint multiplied by a factor covering the cell.
        expand=iv.exp(abs(p.real)*w)*abs(p)**(DEG+1)*w**(DEG+1)/factorial(DEG+1)
        powers=[p**k for k in range(DEG+1)]
        weights=[(p+c)**k for k in range(3)]
        data_rows.append((z,rho,powers,weights,[expand*abs(q) for q in weights]))
    return data_rows,err
P=[I(1)]+[I([0,(w**j).b]) if j%2==0 else I([-(w**j).b,(w**j).b]) for j in range(1,DEG+1)]
P=[p/factorial(j) for j,p in enumerate(P)]
# Precompute exact multipliers to avoid expensive repeated powers.
def build(Q,kind):
    rows,err=setup(Q,kind)
    return [[z,rho,[[p*w for p in powers] for w in weights],errs] for z,rho,powers,weights,errs in rows],err
def jets(rows,err):
    out=[[I(0) for _ in P] for _ in range(3)]
    rem=[I(0)]*3
    for row in rows:
        z,rho,table,est=row
        row[0]=z*rho
        az=abs(z)
        for i in range(3):
            rem[i]+=az*est[i]
            for j,p in enumerate(table[i]):out[i][j]+=(z*p).real
    # Subtract the midpoint logarithmic slope before evaluating the interval
    # polynomials. This is the same algebra as a local exponential tilt,
    # and changes neither logarithmic curvature nor the sign comparison.
    assert out[0][0].a>0
    beta=(out[1][0]/out[0][0]).mid
    a0,a1,a2=out
    rows=[a0,[v-beta*u for u,v in zip(a0,a1)],[v-2*beta*u+beta**2*w for w,u,v in zip(a0,a1,a2)]]
    ae=err.b
    r0,r1,r2=[e+ae for e in rem]
    errors=[r0,r1+abs(beta)*r0,r2+2*abs(beta)*r1+beta**2*r0]
    return [sum(p*b for p,b in zip(P,row))+I([-e.b,e.b]) for row,e in zip(rows,errors)]
X,eX=build(qx,True); T,eT=build(qt,False)
maxdiff=I('-1e30'); minFX=I('100'); minFT=I('100')
for i in range(N):
    x,u,v=jets(X,eX)
    y,p,q=jets(T,eT)
    assert x.a>0 and y.a>0,(i,x,y,eX,eT)
    diff=q/y-(p/y)**2-v/x+(u/x)**2
    assert diff.b<0,(i,diff)
    maxdiff=max(maxdiff,diff.b);minFX=min(minFX,x.a);minFT=min(minFT,y.a)
print(lo,hi,'Q_X,Q_T',qx,qt,'cells',N,'diff upper',maxdiff,'positive unnormalized scaled CDF floors',minFX,minFT,'quadrature intervals',eX,eT)
