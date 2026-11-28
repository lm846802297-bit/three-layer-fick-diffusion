
# Compact-window contour certification with signed first-order parameter response.
from pathlib import Path
from sys import argv
# Root isolation and physical transfer definitions.
oldargv=argv[:]

# Each occurrence of either independent wall rate is enclosed by the same interval.
# Equal-rate endpoint corners bracket eigenvalues by min-max monotonicity.
from mpmath import iv
iv.dps=55
I=iv.mpf
pi=iv.pi
l,n=I(5)/66,I(8)/33
from sys import argv
ends=['3.9','4.1']
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

A0,B0,L0=a,b,la
kap=I(4); ka=kap; ends=['4','4']
aa=root(0,'3.60','3.94')**2
bb=root(1,'4.45','4.72')**2
ll=root(2,'2.23','2.62')**2
from fractions import Fraction
from math import comb, factorial, ceil, sqrt
def sy(v):return I([-v.b,v.b])
def ampdata(q):
    def pair(v):
        w=iv.exp(v); return (w+1/w)/2,(w-1/w)/2
    C,S=pair(l*q); c,s=pair(n*q); U,V=pair(q/2)
    y=C+4*S/(5*q); z=5*S+4*C/q
    Y=c*y+2*s*z; H=s*y/2+c*z
    R=U+4*V/q; J=V+4*U/q
    Y1=(c*S/5+2*s*C)/q; H1=(s*S/10+c*C)/q
    R1=V/q; J1=U/q
    return (R,Y),(0,Y1),(R1,0),Y*J+R*H,Y1*J+R*H1,Y*J1+R1*H,Y1*J1+R1*H1
def nodes(Q,M,j):
    h=I(1)/2
    r0=[bb,aa][j]; r=[B0,A0][j]
    pars=[sy(I('.1')),sy(I('.1')),L0-ll,r-r0]
    rows=[]
    for k in range(M+1):
        q=iv.mpc(Q,k*h); p=q*q
        Y,Yl,Yr,D,Dl,Dr,Dlr=ampdata(q)
        old=(h if k==0 else 2*h)/iv.pi*Y[j]/D*(p+ll)/(p+r0)
        U,V=Yl[j]/Y[j],Yr[j]/Y[j]
        G,H,J=Dl/D,Dr/D,Dlr/D
        S,W=1/(p+ll),1/(p+r0)
        basis=[U-G,V-H,S,-W]
        y=U*pars[0]+V*pars[1]; z=S*pars[2]; w=W*pars[3]
        lin=G*pars[0]+H*pars[1]; quad=J*pars[0]*pars[1]
        d=lin+quad
        resid=(y*z-quad-d*w-(y+z-lin-w)*(d+w+d*w))/((1+d)*(1+w))
        rows.append((p,old,[old*b for b in basis],abs(old)*abs(resid)))
    return rows,pars
def gp(t,Q,k):
    mm=[iv.sqrt(iv.pi/t)/2,1/(2*t)]
    for j in range(2,k+1):mm.append((j-1)*mm[j-2]/(2*t))
    return sum(comb(k,j)*Q**(k-j)*mm[j] for j in range(k+1))
def block(lo,hi,qx,qt,step):
    ncell=(Fraction(hi)-Fraction(lo))/Fraction(step)
    assert ncell.denominator==1
    lo,hi,step=I(lo),I(hi),I(step)
    rad=step/2
    ans=[]
    D=3
    for side,Q in enumerate([qx,qt]):
        d=[l+n,I(1)/2][side]
        strip=3; low=Q-strip; high=Q+strip
        assert low>=4
        # Exact exponential envelope from transfer-matrix leading term.
        B=(1+iv.exp(-2*l*low))*(1+iv.exp(-2*n*low))*(1+iv.exp(-low))
        den=I(33)/16*(2-B)
        assert den>0
        num=[(1+I('4.1')/low)*(1+iv.exp(-low))/2,
             I(11)/4*(1+I('.82')/low)*(1+iv.exp(-2*l*low))*(1+iv.exp(-2*n*low))][side]
        # Re(p)>0 is not required; |p+r|>=|q|^2 since |q +/- i sqrt(r)| product
        # may be smaller, but Re(q)>=low and |(p+l)/(p+r)| <= 1+(r-l)/|p+r|.
        # Min of |q^2+r| over Im(q) is 2*Re(q)*sqrt(r) if r>=Re(q)^2,
        # otherwise Re(q)^2+r. Both are at least low^2.
        C=(num/den*(1+24/I(low)**2))/iv.pi
        assert C<100
        h=I('.5'); M=ceil(2*sqrt(55/float(lo.a))); Y=(M+1)*h
        rows,pars=nodes(Q,M,side)
        ali=iv.exp(-2*iv.pi*strip/h)
        out=[]
        for i in range(int(ncell)):
            t=lo+(i+I('.5'))*step
            beta=(d*d/(4*t*t)-1/(2*t)).mid
            expscale=d*Q-t*Q**2
            coeff=[[I(0)]*(D+1) for _ in range(3)]
            sensitivity=[[[I(0)]*(D+1) for _ in pars] for _ in range(3)]
            er=[I(0)]*3
            w=sy(rad)
            for p,o,wts,res in rows:
                ex=iv.exp(p*t+expscale); v=ex*o
                s=p-beta
                ew=[ex*a for a in wts]
                ab=abs(p); mag=abs(s)
                e=ab*rad
                fact=iv.exp(e)
                rem=fact*e**(D+1)/factorial(D+1)
                bound=abs(v)+sum(abs(a)*abs(b) for a,b in zip(ew,pars))
                # All later operations use centered jets of U, with no derivative of the display scale.
                for j in range(3):
                    er[j]+=mag**j*(bound*rem+abs(ex)*res*fact)
                    power=s**j
                    for u in range(D+1):
                        coeff[j][u]+=(v*power).real
                        for k in range(4):sensitivity[j][k][u]+=(ew[k]*power).real
                        power*=p
            jts=[]
            for j in range(3):
                # |p-beta|^j <= (|q|^2+|beta|)^j on both strip edges.
                tmin=t-rad; tmax=t+rad
                envelope=4*C*iv.exp(expscale+tmax*high**2-d*low)
                aerr=envelope*sum(comb(j,k)*abs(beta)**(j-k)*gp(tmin,high,2*k) for k in range(j+1))*ali/(1-ali)
                rr=((Q+Y+h)/(Q+Y))**(2*j)*iv.exp(-tmin*(2*Y*h+h*h))
                assert rr<1
                terr=2*h*C*iv.exp(expscale+tmax*Q**2-d*Q-tmin*Y**2)*((Q+Y)**2+abs(beta))**j/(1-rr)
                row=[coeff[j][u]+sum(p*sensitivity[j][k][u] for k,p in enumerate(pars)) for u in range(D+1)]
                v=sum(row[u]*w**u/factorial(u) for u in range(D+1))+sy(er[j]+aerr+terr)
                jts.append(v)
            out.append(jts)
        ans.append(out)
    low=I('inf')
    for i,(X,T) in enumerate(zip(*ans)):
        for V in [X,T]:assert V[0]>0,(i,V)
        h=lambda V: V[2]/V[0]-(V[1]/V[0])**2
        v=h(X)-h(T)
        low=min(low,v.a)
    print('block',lo,hi,'lower',low,flush=True)
    return low

if __name__=='__main__':
    if len(oldargv)>1:
        assert block(*oldargv[1:3],*map(int,oldargv[3:5]),oldargv[5])>0
    else:
        for args in [('.035','.062',8,9,'.0005'),('.015','.035',10,12,'.00025'),('.008','.015',22,29,'.0001'),('.005','.008',32,45,'.00005'),('.003','.005',43,66,'.000025'),('.002','.003',63,100,'.00001'),('.0015','.002',90,143,'.000005')]:
            assert block(*args)>0
