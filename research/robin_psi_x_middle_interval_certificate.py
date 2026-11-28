"""70-digit interval certificate for N_X<0 on [0.0015, 0.010].
The exact residues are obtained from the full Robin characteristic determinant.
The infinite tail is bounded analytically; no fitted or truncated-tail extrapolation is used.
"""
import math, platform, sys, time, os, json
import mpmath as mp
from mpmath import iv
import numpy as np
from scipy.optimize import brentq

HERE=os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE,"robin_psi_x_middle_interval_certificate_input.json"),encoding="utf-8") as _f:
    CFG=json.load(_f)
mp.mp.dps=int(CFG["interval_arithmetic"]["mpmath_digits"])
iv.dps=int(CFG["interval_arithmetic"]["interval_digits"])
I=iv.mpf
kap=I(CFG["benchmark"]["kappa"])
w1=I(CFG["benchmark"]["w1"]); w2=I(CFG["benchmark"]["w2"])
wr=I(CFG["benchmark"]["right_half_width"])
D1=I(CFG["benchmark"]["D1"]); D2=I(CFG["benchmark"]["D2"])
DR=I(CFG["benchmark"]["D_right"])
N=int(CFG["spectral"]["finite_residue_indices"][1])
assert int(CFG["spectral"]["global_roots_required"])==N+2
ROOT_HALF_WIDTH=mp.mpf(CFG["spectral"]["root_bracket_half_width"])
ROOT_DIGITS=int(CFG["spectral"]["root_bracket_digits"])


def mat_pair(k,h,D):
    sd=iv.sqrt(D); a=h/sd; z=k*a
    C=iv.cos(z); S=iv.sin(z); Cp=-a*S; Sp=a*C
    M=[[C,S/(sd*k)],[-sd*k*S,C]]
    Mp=[[Cp,(Sp*k-S)/(sd*k*k)],[-sd*(S+k*Sp),Cp]]
    return M,Mp

def mm(A,B):
    return [[A[0][0]*B[0][0]+A[0][1]*B[1][0],A[0][0]*B[0][1]+A[0][1]*B[1][1]],
            [A[1][0]*B[0][0]+A[1][1]*B[1][0],A[1][0]*B[0][1]+A[1][1]*B[1][1]]]

def madd(A,B): return [[A[i][j]+B[i][j] for j in range(2)] for i in range(2)]
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
def mm_mp(A,B): return [[sum(A[i][q]*B[q][j] for q in range(2)) for j in range(2)] for i in range(2)]
def fullF_mp(k):
    M=mm_mp(mat_mp(k,mp.mpf(1)/2,1),mm_mp(mat_mp(k,mp.mpf(4)/33,mp.mpf(1)/4),mat_mp(k,mp.mpf(25)/66,25)))
    y=M[0][0]+4*M[0][1]; p=M[1][0]+4*M[1][1]
    return p+4*y
def mat_float(k,h,D):
    z=k*h/math.sqrt(D); c=math.cos(z); s=math.sin(z)
    return np.array([[c,s/(math.sqrt(D)*k)],[-math.sqrt(D)*k*s,c]])
def fullF_float(k):
    M=mat_float(k,.5,1)@mat_float(k,4/33,.25)@mat_float(k,25/66,25)
    y=M[0,0]+4*M[0,1]; p=M[1,0]+4*M[1,1]
    return p+4*y

def scan(fun,kmax,step,count):
    grid=np.arange(1e-5,kmax+step,step); vals=np.array([fun(float(x)) for x in grid]); out=[]
    for a,b,fa,fb in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]):
        if fa*fb<0:
            x=brentq(fun,a,b,xtol=1e-13,rtol=1e-14)
            if not out or x-out[-1]>1e-7: out.append(x)
            if len(out)==count: break
    if len(out)!=count: raise RuntimeError((len(out),count))
    return out

def node_count_interval(kI):
    y=I(1); p=I(4); total=0
    specs=[(I(25)/66,I(25)),(I(4)/33,I(1)/4),(I(1)/2,I(1))]
    for h,D in specs:
        sd=iv.sqrt(D)
        u=p/(sd*kI)
        phi=iv.atan2(u,y)
        theta=kI*h/sd
        # Midpoint only selects candidate integer labels; interval tests certify them.
        km=(mp.mpf(kI.a)+mp.mpf(kI.b))/2
        ym=mp.mpf(y.mid); pm=mp.mpf(p.mid); Dm=mp.mpf(D.mid); hm=mp.mpf(h.mid)
        phim=mp.atan2(pm/(mp.sqrt(Dm)*km),ym)
        thetam=km*hm/mp.sqrt(Dm)
        low=(phim-mp.pi/2-thetam)/mp.pi
        high=(phim-mp.pi/2)/mp.pi
        m0=int(mp.floor(low))-1; m1=int(mp.ceil(high))+1
        for m in range(m0,m1+1):
            z=phi-iv.pi/2-m*iv.pi
            if z.a>0 and z.b<theta.a:
                total += 1
            elif z.b<0 or z.a>theta.b:
                pass
            else:
                raise RuntimeError(('ambiguous Sturm zero',kI,h,D,phi,theta,z,m))
        M=mat_pair(kI,h,D)[0]
        y,p=M[0][0]*y+M[0][1]*p,M[1][0]*y+M[1][1]*p
    return total

def isolate(fun_iv,root,digits=ROOT_DIGITS):
    rad=ROOT_HALF_WIDTH
    lo=mp.nstr(root-rad,digits+12); hi=mp.nstr(root+rad,digits+12)
    flo=fun_iv(I(lo)); fhi=fun_iv(I(hi))
    assert (flo.b<0 and fhi.a>0) or (flo.a>0 and fhi.b<0),(lo,hi,flo,fhi)
    return I([lo,hi])

start=time.time()
kg_seed=scan(fullF_float,float(CFG["spectral"]["root_scan_kmax"]),float(CFG["spectral"]["root_scan_step"]),N+2)
kg=[mp.findroot(fullF_mp,(mp.mpf(x)-mp.mpf('1e-5'),mp.mpf(x)+mp.mpf('1e-5'))) for x in kg_seed]
kI=[isolate(lambda x:full_pair_iv(x)[0],x) for x in kg]
nodes=[node_count_interval(ki) for ki in kI]
assert nodes==list(range(N+2)), nodes
assert all((full_pair_iv(I(mp.nstr(x-ROOT_HALF_WIDTH,70)))[0].b<0 and full_pair_iv(I(mp.nstr(x+ROOT_HALF_WIDTH,70)))[0].a>0) or (full_pair_iv(I(mp.nstr(x-ROOT_HALF_WIDTH,70)))[0].a>0 and full_pair_iv(I(mp.nstr(x+ROOT_HALF_WIDTH,70)))[0].b<0) for x in kg)
lam=[x*x for x in kI]; lam1=lam[0]; c=lam[1]

def PGp(k):
    _,dp=full_pair_iv(k)
    return -dp/(2*k*24)
def PR(k): return right_pair_iv(k)[0]/3

# R1 is the right-half first positive root, reused from the earlier Q30 interval certificate.
kR1=I(CFG["spectral"]["right_first_root_wave_number_bracket"])
R1=kR1*kR1
PGpc=PGp(kI[1]); PRc=PR(kI[1])
C=(1-c/R1)/PRc*(c*PGpc)/(1-c/lam1)
assert abs(C).b<1
eta=[]; B=[]
for i in range(2,N+1): # B_i for i=2..N, global lambda_{i+1}; kI index i
    ln=lam[i]
    ei=ln-c
    bi=C*PR(kI[i])/(1-ln/R1)*(1-ln/lam1)*(1-ln/c)/PGp(kI[i])
    eta.append(ei); B.append(bi)
assert B[0].a>0 and B[1].a>0
for i in range(N-2): assert eta[i+1].a>eta[i].b

# Tail bounds from Neumann bracketing and the Lagrange residue identity.
# Neumann bracketing gives lambda_n >= pi^2 (n-3)^2/(9/11)^2.
# For i>=81, n=i+1, this implies eta_i=lambda_{i+1}-c >= 14 i^2.
coef=iv.pi**2*I(121)/81
assert coef.a>I('14.74')
g81=I('14.74')*I(79)**2-I('15.81')-I(14)*I(81)**2
stepg=I('2')*I('14.74')*I(81)-3*I('14.74')-I(14)*I(163)
assert g81.a>0 and stepg.a>0
# Dirichlet trial subspaces give lambda_{i+1} <=25*pi^2*(i+1)^2 <=100*pi^2*i^2.
# For k>=10, Lagrange identity and first-layer norm yield |P_G'(-lambda)|>0.029/(24*36).
# The residue formula then gives |B_i|<=3e4*lambda_{i+1}<=3e7*i^2.
PGp_lower=I('0.029')/(24*36)
residue_coeff=(24*36/I('0.029'))*(I('1.4')/3)*(2*R1.b/(lam1.a*c.a))
assert residue_coeff.b<I('3e4')
tau=I(CFG["time_certification"]["tau"]); Cb=I(CFG["infinite_tail_bounds"]["residue_constant"]); Ce=I(100)*iv.pi**2; al=I(CFG["infinite_tail_bounds"]["rate_lower_constant"])
def sum_tail(power):
    j=N+1
    term=I(j)**power*iv.exp(-al*tau*j*j)
    ratio=(I(j+1)/j)**power*iv.exp(-al*tau*(2*j+1))
    assert ratio.b<1
    return term/(1-ratio)
tails={m:Cb*Ce**m*sum_tail(2*m+2) for m in range(3)}
Ftail=(Cb/al)*sum_tail(0)

M=int(CFG["interval_arithmetic"]["taylor_degree"])
max_m=M+3
Sbounds=[]
for m in range(max_m+1):
    finite=I(0)
    for i in range(N-1):
        finite += abs(B[i])*eta[i]**m*iv.exp(-eta[i]*tau)
    tail=Cb*Ce**m*sum_tail(2*m+2)
    Sbounds.append(finite.b+tail.b)

def Fbound(m):
    return I(1) if m==0 else I(Sbounds[m-1])

def Nderivative_bound(n):
    total=I(0)
    for k in range(n+1):
        total += math.comb(n,k)*I(Sbounds[k+2])*Fbound(n-k)
        total += math.comb(n,k)*I(Sbounds[k+1])*I(Sbounds[n-k])
    return total.b
Nrem_bound=Nderivative_bound(M+1)

def taylor_box(tlo,thi):
    center=(tlo+thi)/2
    radius=(thi-tlo)/2
    center_err=mp.mpf(CFG["interval_arithmetic"]["time_center_enclosure_radius"])
    tc=I([mp.nstr(center-center_err,75),mp.nstr(center+center_err,75)])
    fder=[I(0) for _ in range(M+4)]
    F0=I(1)
    for i in range(N-1):
        e=eta[i]; b=B[i]; ex=iv.exp(-e*tc)
        F0 -= b/e*ex
        powe=I(1)
        for k in range(M+4):
            if k>0: powe *= -e
            fder[k] += b*powe*ex
    for k in range(M+4):
        T=Cb*Ce**k*sum_tail(2*k+2)
        fder[k] += I([-T.b,T.b])
    F0 += I([-Ftail.b,Ftail.b])
    Fder=[F0]+[fder[m-1] for m in range(1,M+2)]
    coeff=[]
    for n in range(M+1):
        an=I(0); bn=I(0)
        for k in range(n+1):
            an += math.comb(n,k)*fder[k+2]*Fder[n-k]
            bn += math.comb(n,k)*fder[k+1]*fder[n-k]
        coeff.append((an-bn)/mp.factorial(n))
    h=I([-(radius+center_err),radius+center_err])
    poly=coeff[M]
    for n in range(M-1,-1,-1):
        poly=poly*h+coeff[n]
    rem=I(Nrem_bound)*I(mp.nstr(radius+center_err,75))**(M+1)/mp.factorial(M+1)
    return poly+I([-rem.b,rem.b]), rem

def certify(step):
    a=mp.mpf(CFG["time_certification"]["start"]); b=mp.mpf(CFG["time_certification"]["end"]); n=int(mp.ceil((b-a)/step))
    max_upper=None; min_lower=None; max_rem=None
    for j in range(n):
        lo=a+(b-a)*j/n; hi=a+(b-a)*(j+1)/n
        NX,rem=taylor_box(lo,hi)
        if NX.b>=0: raise RuntimeError(('box fail',j,lo,hi,NX,rem))
        max_upper=NX.b if max_upper is None or NX.b>max_upper else max_upper
        min_lower=NX.a if min_lower is None or NX.a<min_lower else min_lower
        max_rem=rem.b if max_rem is None or rem.b>max_rem else max_rem
    return n,max_upper,min_lower,max_rem

print('N',N,'root max',mp.nstr(kg[-1],15),'tail bounds',tails,'Ftail',Ftail)
print('S_m bounds m=0..11',[mp.nstr(x,12) for x in Sbounds])
print('Taylor degree',M,'uniform N derivative remainder bound',mp.nstr(Nrem_bound,15))
for step in [CFG["time_certification"]["maximum_box_width"], '0.00001','0.000005']:
    try:
        n,upper,lower,maxrem=certify(mp.mpf(step))
        assert upper < I(CFG["acceptance"]["N_X_upper_bound"]).a
        assert n==int(CFG["time_certification"]["expected_boxes"])
        print('PASS step',step,'boxes',n,'sup NX upper',upper,'inf NX lower',lower,'max Taylor remainder',maxrem)
        break
    except RuntimeError as e:
        print('FAIL step',step,str(e)[:500])
else: raise RuntimeError('all box partitions failed')
print('root lambda1,c,R1',lam1,c,R1)
print('C',C,'B2',B[0],'B3',B[1],'eta2',eta[0],'eta3',eta[1])
print('global interval Sturm node counts',len(nodes),'verified; counts 0..81')
print('Neumann lower coefficient pi^2*(11/9)^2',coef)
print('eta_i >= 14*i^2 proof checks at i=81 and monotone difference',g81,stepg)
print('residue coefficient upper for |B_i|<=C*lambda',residue_coeff)
print('S_m bounds m=0..11',[mp.nstr(x,12) for x in Sbounds])
print('Taylor degree',M,'uniform N derivative remainder bound',mp.nstr(Nrem_bound,15))
print('time boxes',n,'max box width',mp.nstr((mp.mpf(CFG["time_certification"]["end"])-mp.mpf(CFG["time_certification"]["start"]))/n,12))
print('runtime_seconds',time.time()-start,'platform',platform.platform(),'python',sys.version.split()[0],'mpmath',mp.__version__,'numpy',np.__version__,'scipy',__import__('scipy').__version__,'cpu_count',os.cpu_count(),'gpu=not_used')
