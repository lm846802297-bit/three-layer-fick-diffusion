# Each occurrence of either independent wall rate is enclosed by the same interval.
# Equal-rate endpoint corners bracket eigenvalues by min-max monotonicity.
from mpmath import iv
from sys import argv
iv.dps=55
I=iv.mpf
from sys import argv
ends=argv[1:] or ['3.99999','4.00001']
assert len(ends)==2
kap=I(ends)
ka=kap
def data(q,ctx):
    S=lambda x:(ctx.exp(x)+ctx.exp(-x))/2
    H=lambda x:(ctx.exp(x)-ctx.exp(-x))/2
    l=ctx.mpf(5)/66; n=ctx.mpf(8)/33
    c,s,C,t,U,V=S(l*q),H(l*q),S(n*q),H(n*q),S(q/2),H(q/2)
    y=c+ka/5/q*s; z=5*s+ka/q*c
    L=C*y+2*t*z; BL=t*y/2+C*z
    R=U+ka/q*V; BR=V+ka/q*U
    return L,R,L*BR+R*BL
def spec(x):return data(iv.mpc(0,x),iv)[2].imag
def rs(x):return data(iv.mpc(0,x),iv)[1].real
def ls(x):return data(iv.mpc(0,x),iv)[0].real
def root0(fun,x,y):
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

def root(fun,x,y):
    global ka
    ka=kap
    assert (fun(I(x))*fun(I(y))).b<0
    ka=I(ends[0]); z0=root0(fun,x,y)
    ka=I(ends[1]); z1=root0(fun,x,y)
    ka=kap
    return I([z0.a,z1.b])
assert 0<kap<I('4.01')
a=root(ls,'3.77','3.78'); b=root(rs,'4.5','4.6')
la=root(spec,'2.41','2.43'); c=root(spec,'3.9','4')
assert la<a<c<b<21



from math import comb,factorial,ceil,sqrt
from fractions import Fraction as Fr
rho=root(rs,'10.17','10.18')-c
l,n=I(5)/66,I(8)/33
ep=I('.035'); t0=I('.0015')
def sym(x):return I([-x.b,x.b])
def mom(t,N):
    out=[iv.sqrt(iv.pi/t)/2,1/(2*t)]
    for j in range(2,N+1): out.append((j-1)*out[j-2]/(2*t))
    return out
def gaussian_poly(t,Q,N):
    mm=mom(t,N)
    return sum(comb(N,j)*Q**(N-j)*mm[j] for j in range(N+1))
def block(lo,hi,Q,step,jets,kind='X',score=0,center=False):
    # Holomorphic strip in y, and absolute Gaussian bounds.
    tlo,thi,dd=map(I,(lo,hi,step))
    rad=dd/2
    A=2 if kind=='Z' else 3; h=I(1)/4; low=Q-A; high=Q+A
    d=(l+n if kind=='X' else I(1)/2) if kind!='Z' else I(2)/11
    assert low>= (1 if kind=='Z' else 4)
    bb=(1+iv.exp(-2*l*low))*(1+iv.exp(-2*n*low))*(1+iv.exp(-low))
    if kind!='Z': assert I(33)/16*(2-bb)>(11+I('8.83')/low)/100
    assert 2*(11+I('8.83')/low)/(1-iv.exp(-low))<100
    # This bound includes 1/pi, using pi>1.
    k=100; fac=1+I(21)/low**2
    M=ceil(4*sqrt(110/float(lo))); Y=(M+1)*h
    aa=iv.exp(-2*iv.pi*A/h)
    errs={}
    for j in jets:
        m=2*j+1+(2 if score else 0)
        ext=1+I(116)/low**2 if score else 1
        ali=4*k*ext*fac**(j+1)*iv.exp(thi*(high**2+c)-d*low)*gaussian_poly(tlo,high,m)*aa/(1-aa)
        rr=((Q+Y+h)/(Q+Y))**m*iv.exp(-tlo*(2*Y*h+h*h))
        assert rr<1
        tail=2*h*k*ext*fac**(j+1)*iv.exp(thi*(Q**2+c)-d*Q-tlo*Y**2)*(Q+Y)**m/(1-rr)
        errs[j]=ali+tail
    D=5
    rows=[]
    for k in range(M+1):
        q=iv.mpc(Q,k*h); p=q*q; L,R,B=data(q,iv)
        z=(R/B*(p+la)/(p+b)) if kind=='X' else ((L/B*(p+la)/(p+a)) if kind=='T' else (q*L/R*(p+b)/(p+a)))
        s=p+c
        v=(1 if k==0 else 2)*h/iv.pi*z*(s+score if score else 1)
        ab=abs(s)
        env=[ab**j*(ab*rad)**(D+1)*iv.exp(ab*rad)/factorial(D+1) for j in jets]
        rows.append([v,p,[s**j for j in range(max(jets)+D+1)],env])
    count=1 if Fr(step)==0 else (Fr(hi)-Fr(lo))/Fr(step)
    assert isinstance(count,int) or count.denominator==1
    for i in range(int(count)):
        coeff=[I(0)]*(max(jets)+D+1)
        rem=[I(0) for j in jets]
        tc=tlo+dd*i+rad
        ec=iv.exp(c*tc)
        for v,p,pows,env in rows:
            v=v*iv.exp(p*tc)*ec
            for j,p in enumerate(pows):
                coeff[j]+=(v*p).real
            av=abs(v)
            for k,weight in enumerate(env):
                rem[k]+=av*weight
        w=sym(rad)
        if center:
            assert jets==(0,1,2)
            beta=d**2/(4*tc**2)-1/(2*tc)
            beta=beta.mid
            errs2=[rem[v]+errs[v] for v in jets]
            def jetc(j):
                combs=[comb(j,i)*(-beta)**(j-i) for i in range(j+1)]
                tay=[sum(combs[i]*coeff[i+k] for i in range(j+1)) for k in range(D+1)]
                rr=sum(abs(combs[i])*errs2[i] for i in range(j+1))
                return sum(tay[k]*w**k/factorial(k) for k in range(D+1))+sym(rr)
            yield [jetc(j) for j in jets]
        else:
            yield [sum(coeff[j+k]*w**k/factorial(k) for k in range(D+1))+sym(rem[i]+errs[j]) for i,j in enumerate(jets)]
def one(t,Q,j,kind='X'):
    return next(block(t,t,Q,'0',j,kind))
F,F2=one('.07',7,(0,2))
m=F2/F
assert m<0
rho=I(85); Mz=I(95)
print('m M rho',m,Mz,rho,flush=True)
# Short omitted slice: saddle bound |F^(j)| with |e^(dq)Y/D*(p+la)/(p+b)|<9.
# Each term C t^(-power) exp(ct-A/t) increases up to t0 since A/t0>power.
d=l+n
Q=d/(2*t0); AA=d*d/4
assert AA/t0>I('4.5')
gh=iv.sqrt(iv.pi)/2
fac=18/iv.pi*iv.exp(c*t0-AA/t0)/iv.sqrt(t0)*gh
g0=fac
u=Q**2+21
g2=fac*(u*u+u/t0+I(3)/(4*t0*t0))
early=t0*iv.exp(Mz*t0)*(g2+abs(m)*g0)
overall=-early.b
print('early absolute bound',early,flush=True)
for lo,hi,Q,step in [('.0015','.005',48,'.0001'),('.005','.015',22,'.0002'),('.015','.035',9,'.00025')]:
    dd=I(step); subtotal=I(0)
    for i,(F,F2) in enumerate(block(lo,hi,Q,step,(0,2))):
        v=I(lo)+I([i,i+1])*dd
        B=m*F-F2
        val=iv.exp(rho*v)*I([max(I(0),B.a),max(I(0),B.b)])-iv.exp(Mz*v)*I([max(I(0),-B.b),max(I(0),-B.a)])
        subtotal+=dd*val
    overall+=subtotal.a
    print('integral block',lo,hi,subtotal,flush=True)
print('pivot lower',overall,flush=True)
assert overall>0
for lo,hi,Q,step in [('.035','.07',7,'.0005'),('.07','.10',7,'.0005')]:
    upper=I('-inf'); lower=I('inf')
    for fp,fpp in block(lo,hi,Q,step,(2,3)):
        upper=max(upper,fp.b);lower=min(lower,fpp.a)
    print('suffix bounds',lo,hi,upper,lower,flush=True)
    assert upper<0<lower

# Direct compact-window score bounds, without assuming a curvature theorem.
slow=max(j[0].b for j in block('.035','.07',3,'.0002',(0,),'Z',85))
fast=min(j[0].a for j in block('.035','.07',3,'.0002',(0,),'Z',95))
print('Z compact bound zprime+85z upper, zprime+95z lower:',slow,fast,flush=True)
assert slow<0<fast

# Join the pivot conclusion to the parameter-uniform analytic short-time slice.
for lo,hi,qx,qt,step in [('.0015','.002',91,143,'.000005'),('.002','.003',65,102,'.000005'),('.003','.005',41,64,'.00001'),('.005','.015',22,32,'.00005'),('.015','.035',9,12,'.0002'),('.035','.07',7,8,'.00025')]:
    maxdiff=I('-inf')
    for xx,tt in zip(block(lo,hi,qx,step,(0,1,2),'X',center=True),block(lo,hi,qt,step,(0,1,2),'T',center=True)):
        assert min(xx[0].a,tt[0].a)>0,(lo,hi,xx[0],tt[0])
        def hh(j): return j[2]/j[0]-(j[1]/j[0])**2
        maxdiff=max(maxdiff,(hh(tt)-hh(xx)).b)
    print('compact hTprime-hXprime upper',lo,hi,maxdiff,flush=True)
    assert maxdiff<0
