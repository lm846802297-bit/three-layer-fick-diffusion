import math, platform, sys, time
import mpmath as mp
from mpmath import iv
import numpy as np
from scipy.optimize import brentq

mp.mp.dps=90
iv.dps=70
I=iv.mpf
kap=I(4)
w1,w2,wr=I(25)/66,I(4)/33,I(1)/2
D1,D2,DR=I(25),I(1)/4,I(1)
N=100

def mat_pair(k,h,D):
    sd=iv.sqrt(D); a=h/sd; z=k*a
    C=iv.cos(z); S=iv.sin(z)
    Cp=-a*S; Sp=a*C
    M=[[C,S/(sd*k)],[-sd*k*S,C]]
    Mp=[[Cp,(Sp*k-S)/(sd*k*k)],[-sd*(S+k*Sp),Cp]]
    return M,Mp

def mm(A,B):
    return [[A[0][0]*B[0][0]+A[0][1]*B[1][0],A[0][0]*B[0][1]+A[0][1]*B[1][1]],
            [A[1][0]*B[0][0]+A[1][1]*B[1][0],A[1][0]*B[0][1]+A[1][1]*B[1][1]]]

def madd(A,B):
    return [[A[i][j]+B[i][j] for j in range(2)] for i in range(2)]

def pairmul(P,Q):
    A,Ap=P; B,Bp=Q
    return mm(A,B),madd(mm(Ap,B),mm(A,Bp))

def full_pair_iv(k):
    MR=mat_pair(k,wr,DR); M2=mat_pair(k,w2,D2); M1=mat_pair(k,w1,D1)
    A,Ap=pairmul(pairmul(MR,M2),M1)
    y=A[0][0]+A[0][1]*kap; p=A[1][0]+A[1][1]*kap
    yp=Ap[0][0]+Ap[0][1]*kap; pp=Ap[1][0]+Ap[1][1]*kap
    return p+kap*y, pp+kap*yp

def right_pair_iv(k):
    Y=iv.cos(k/2)+kap/k*iv.sin(k/2)
    Yp=-iv.sin(k/2)/2+kap*(iv.cos(k/2)/(2*k)-iv.sin(k/2)/(k*k))
    return Y,Yp

def mat_mp(k,h,D):
    z=k*h/mp.sqrt(D); c=mp.cos(z); s=mp.sin(z)
    return [[c,s/(mp.sqrt(D)*k)],[-mp.sqrt(D)*k*s,c]]

def mm_mp(A,B):
    return [[sum(A[i][q]*B[q][j] for q in range(2)) for j in range(2)] for i in range(2)]

def fullF_mp(k):
    M=mm_mp(mat_mp(k,mp.mpf(1)/2,1),mm_mp(mat_mp(k,mp.mpf(4)/33,mp.mpf(1)/4),mat_mp(k,mp.mpf(25)/66,25)))
    y=M[0][0]+4*M[0][1]; p=M[1][0]+4*M[1][1]
    return p+4*y

def rightY_mp(k):
    return mp.cos(k/2)+4/k*mp.sin(k/2)

def mat_float(k,h,D):
    z=k*h/math.sqrt(D); c=math.cos(z); s=math.sin(z)
    return np.array([[c,s/(math.sqrt(D)*k)],[-math.sqrt(D)*k*s,c]])

def fullF_float(k):
    M=mat_float(k,.5,1)@mat_float(k,4/33,.25)@mat_float(k,25/66,25)
    y=M[0,0]+4*M[0,1]; p=M[1,0]+4*M[1,1]
    return p+4*y

def rightY_float(k):
    return math.cos(k/2)+4/k*math.sin(k/2)

def node_count(k):
    y=mp.mpf(1); p=mp.mpf(4); count=0
    for h,D in [(mp.mpf(25)/66,mp.mpf(25)),(mp.mpf(4)/33,mp.mpf(1)/4),(mp.mpf(1)/2,mp.mpf(1))]:
        phi=mp.atan2(p/(mp.sqrt(D)*k),y)
        theta=k*h/mp.sqrt(D)
        base=phi+mp.pi/2
        mmin=mp.floor(-base/mp.pi)+1
        mmax=mp.ceil((theta-base)/mp.pi)-1
        count += max(0,int(mmax-mmin+1))
        M=mat_mp(k,h,D)
        y,p=M[0][0]*y+M[0][1]*p,M[1][0]*y+M[1][1]*p
    return count

def scan(fun,kmax,step,count):
    grid=np.arange(1e-5,kmax+step,step); vals=np.array([fun(float(x)) for x in grid]); out=[]
    for a,b,fa,fb in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]):
        if fa*fb<0:
            x=brentq(fun,a,b,xtol=1e-13,rtol=1e-14)
            if not out or x-out[-1]>1e-7: out.append(x)
            if len(out)==count: break
    if len(out)!=count: raise RuntimeError((len(out),count))
    return out

def isolate(fun_iv,root,digits=35):
    rad=mp.mpf(10)**(-digits); lo=mp.nstr(root-rad,digits+12); hi=mp.nstr(root+rad,digits+12)
    flo=fun_iv(I(lo)); fhi=fun_iv(I(hi))
    assert (flo.b<0 and fhi.a>0) or (flo.a>0 and fhi.b<0), (lo,hi,flo,fhi)
    return I([lo,hi])

start=time.time()
kg_seed=scan(fullF_float,500,.05,N+2)
kr_seed=[brentq(rightY_float,(2*j-1)*math.pi+1e-9,2*j*math.pi-1e-9,xtol=1e-13,rtol=1e-14) for j in range(1,N+2)]
kg=[mp.findroot(fullF_mp,(mp.mpf(x)-mp.mpf('1e-5'),mp.mpf(x)+mp.mpf('1e-5'))) for x in kg_seed]
node_counts=[node_count(k) for k in kg]
assert node_counts==list(range(len(kg)))
kr=[mp.findroot(rightY_mp,(mp.mpf(x)-mp.mpf('1e-5'),mp.mpf(x)+mp.mpf('1e-5'))) for x in kr_seed]
lam=[x*x for x in kg]; RR=[x*x for x in kr]
lamI=[isolate(lambda k:full_pair_iv(k)[0],x) for x in kg]
rightKI=[isolate(lambda k:right_pair_iv(k)[0],x) for x in kr]
for j,kI in enumerate(rightKI,start=1):
    assert kI.a>(2*j-1)*iv.pi and kI.b<2*j*iv.pi

c=lamI[1]**2; lam1=lamI[0]**2; R1=rightKI[0]**2

def PGp(kI):
    _,Fp=full_pair_iv(kI)
    return -Fp/(2*kI*24)

def PR(kI):
    return right_pair_iv(kI)[0]/3

PGpc=PGp(lamI[1]); PRc=PR(lamI[1])
C=(1-c/R1)/PRc*(c*PGpc)/(1-c/lam1)
assert abs(C).b<1
eta=[]; B=[]
for i in range(2,N+2):
    ln=lamI[i]**2
    ei=ln-c
    Bi=C*PR(lamI[i])/(1-ln/R1)*(1-ln/lam1)*(1-ln/c)/PGp(lamI[i])
    eta.append(ei); B.append(Bi)
assert B[0].a>0 and B[1].a>0
for i in range(N-1):
    assert eta[i+1].a>eta[i].b

CB=I(3e7); CE=100*iv.pi**2; a=I(2.4)
def tail_sum(power,tmin):
    j=N+1
    term=I(j)**power*iv.exp(-a*tmin*j*j)
    ratio=(I(j+1)/j)**power*iv.exp(-a*tmin*(2*j+1))
    assert ratio.b<1
    return term/(1-ratio)

def tail_f(m,tmin):
    return CB*CE**m*tail_sum(2*m+2,tmin)

def tail_F(tmin):
    return CB/a*tail_sum(0,tmin)

def values(tI,tmin):
    f=I(0); F=I(1); fp=I(0); fpp=I(0); fppp=I(0)
    for i in range(N):
        e=eta[i]; b=B[i]; ex=iv.exp(-e*tI)
        f += b*ex
        F -= b/e*ex
        fp -= e*b*ex
        fpp += e*e*b*ex
        fppp -= e**3*b*ex
    f += I([-tail_f(0,tmin).b,tail_f(0,tmin).b])
    F += I([-tail_F(tmin).b,tail_F(tmin).b])
    fp += I([-tail_f(1,tmin).b,tail_f(1,tmin).b])
    fpp += I([-tail_f(2,tmin).b,tail_f(2,tmin).b])
    fppp += I([-tail_f(3,tmin).b,tail_f(3,tmin).b])
    return f,F,fp,fpp,fppp

def N_and_Np(tI,tmin):
    f,F,fp,fpp,fppp=values(tI,tmin)
    return fpp*F-fp*f, fppp*F-fp*fp

def certify_sign(a0,b0,step,want):
    a0=mp.mpf(a0); b0=mp.mpf(b0); n=int(mp.ceil((b0-a0)/mp.mpf(step)))
    bound=None
    for j in range(n):
        lo=a0+(b0-a0)*j/n; hi=a0+(b0-a0)*(j+1)/n
        tI=I([mp.nstr(lo,45),mp.nstr(hi,45)])
        Nint,_=N_and_Np(tI,I(mp.nstr(lo,45)))
        if want=="negative":
            if Nint.b>=0: raise RuntimeError(("negative fail",lo,hi,Nint))
            if bound is None or Nint.b>bound: bound=Nint.b
        else:
            if Nint.a<=0: raise RuntimeError(("positive fail",lo,hi,Nint))
            if bound is None or Nint.a<bound: bound=Nint.a
    return bound

neg_upper=certify_sign("0.010","0.0144","0.000002","negative")
pos_lower=certify_sign("0.0145","0.020","0.000002","positive")
np_lower=None
for j in range(100):
    lo=mp.mpf('0.0144')+mp.mpf('0.0001')*j/100
    hi=mp.mpf('0.0144')+mp.mpf('0.0001')*(j+1)/100
    tI=I([mp.nstr(lo,45),mp.nstr(hi,45)])
    _,Np=N_and_Np(tI,I(mp.nstr(lo,45)))
    if Np.a<=0: raise RuntimeError(("Nprime fail",lo,hi,Np))
    if np_lower is None or Np.a<np_lower: np_lower=Np.a

t0=I("0.020")
first_fprime=-(eta[0]*B[0]*iv.exp(-eta[0]*t0)+eta[1]*B[1]*iv.exp(-eta[1]*t0))
first_fsecond=eta[0]**2*B[0]*iv.exp(-eta[0]*t0)+eta[1]**2*B[1]*iv.exp(-eta[1]*t0)
tail_fprime=I(0); tail_fsecond=I(0)
for i in range(2,N):
    ex=iv.exp(-eta[i]*t0)
    tail_fprime += max(abs(B[i].a),abs(B[i].b))*max(abs(eta[i].a),abs(eta[i].b))*ex
    tail_fsecond += max(abs(B[i].a),abs(B[i].b))*max(abs(eta[i].a),abs(eta[i].b))**2*ex
tail_fprime += tail_f(1,t0).b
tail_fsecond += tail_f(2,t0).b
assert first_fprime.b<0 and -first_fprime.a>tail_fprime.b
assert first_fsecond.a>tail_fsecond.b

print("ROBIN_PSI_X_LATE_SINGLE_CROSSING_CERTIFICATE: PASS")
print("N modes",N,"interval digits",iv.dps)
print("global Sturm node counts verified",len(node_counts),"roots: n-th root has n-1 interior zeros")
print("right-half roots enclosed in their Sturm brackets ((2j-1)pi,2jpi)")
print("lambda2",c)
print("B2",B[0]); print("B3",B[1])
print("N_X upper on [0.010,0.0144]",neg_upper)
print("N_X lower on [0.0145,0.020]",pos_lower)
print("N_X_prime lower on [0.0144,0.0145]",np_lower)
print("fprime first-pair",first_fprime,"tail abs",tail_fprime)
print("fsecond first-pair",first_fsecond,"tail abs",tail_fsecond)
print("fprime <0 and fsecond >0 for t>=0.020 by exponential ratio monotonicity")
print("runtime",time.time()-start,platform.platform(),sys.version.split()[0],mp.__version__)
sys.exit(0)
