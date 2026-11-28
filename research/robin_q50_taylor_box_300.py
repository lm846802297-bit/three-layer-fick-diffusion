import math, os, platform, sys, time
import mpmath as mp
from mpmath import iv

mp.mp.dps = 90
iv.dps = 70
I = iv.mpf
kap = I(4)
w1, w2, wr = I(25)/66, I(4)/33, I(1)/2
D1, D2, DR = I(25), I(1)/4, I(1)
N = 50

# Interval transfer matrix for state (y, D y').
def mat_iv(k, h, D):
    z = k*h/iv.sqrt(D)
    return [[iv.cos(z), iv.sin(z)/(iv.sqrt(D)*k)],
            [-iv.sqrt(D)*k*iv.sin(z), iv.cos(z)]]
def matmul(A,B):
    return [[A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]],
            [A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]]]
def leftY_iv(k):
    M=matmul(mat_iv(k,w2,D2),mat_iv(k,w1,D1))
    return M[0][0]+M[0][1]*kap
def left_phase_iv(k):
    theta=k*I(5)/66-iv.atan2(I(4)/5,k)
    return iv.atan2(10*iv.sin(theta),iv.cos(theta))+2*k*w2
def rightY_iv(k):
    return iv.cos(k*wr)+kap/k*iv.sin(k*wr)
def rightYp_iv(k):
    return -wr*iv.sin(k*wr)+kap*(wr*iv.cos(k*wr)/k-iv.sin(k*wr)/(k*k))
def fullF_iv(k):
    M=matmul(mat_iv(k,wr,DR),matmul(mat_iv(k,w2,D2),mat_iv(k,w1,D1)))
    y=M[0][0]+M[0][1]*kap
    p=M[1][0]+M[1][1]*kap
    return p+kap*y

# High-precision scalar functions propose narrow root brackets.
def mat_mp(k,h,D):
    z=k*h/mp.sqrt(D); c=mp.cos(z); s=mp.sin(z)
    return [[c,s/(mp.sqrt(D)*k)],[-mp.sqrt(D)*k*s,c]]
def mm(A,B):
    return [[sum(A[i][q]*B[q][j] for q in range(2)) for j in range(2)] for i in range(2)]
def leftY_mp(k):
    M=mm(mat_mp(k,mp.mpf(4)/33,mp.mpf(1)/4),mat_mp(k,mp.mpf(25)/66,mp.mpf(25)))
    return M[0][0]+M[0][1]*4
def rightY_mp(k):
    return mp.cos(k/2)+4/k*mp.sin(k/2)
def fullF_mp(k):
    M=mm(mat_mp(k,mp.mpf(1)/2,1),mm(mat_mp(k,mp.mpf(4)/33,mp.mpf(1)/4),mat_mp(k,mp.mpf(25)/66,25)))
    y=M[0][0]+4*M[0][1]; p=M[1][0]+4*M[1][1]
    return p+4*y

def bracket_root(root, digits=25):
    rad=mp.mpf(10)**(-digits)
    return mp.nstr(root-rad, digits+12), mp.nstr(root+rad, digits+12)
def sign_change(fun,lo,hi,label):
    fl=fun(I(lo)); fh=fun(I(hi))
    assert (fl.b < 0 and fh.a > 0) or (fl.a > 0 and fh.b < 0), (label,fl,fh)
    return I([lo,hi])

start=time.time()
kL1_mp=mp.findroot(leftY_mp, (3.7,3.8))
kR1_mp=mp.findroot(rightY_mp, (4.5,4.7))
kG2_mp=mp.findroot(fullF_mp, (3.9,4.1))
L1k=sign_change(leftY_iv,*bracket_root(kL1_mp),"L1")
phase_lo=left_phase_iv(I(bracket_root(kL1_mp)[0])); phase_hi=left_phase_iv(I(bracket_root(kL1_mp)[1]))
assert phase_lo.b < iv.pi/2 and phase_hi.a > iv.pi/2
R1k=sign_change(rightY_iv,*bracket_root(kR1_mp),"R1")
c_k=sign_change(fullF_iv,*bracket_root(kG2_mp),"global lambda2")
L1=L1k**2; R1=R1k**2; c=c_k**2
assert L1.a < c.a and c.b < R1.a
assert c.b < 5*iv.pi**2

right_ks={1:R1k}
for j in range(2,N+1):
    guess=(2*j-1)*mp.pi + mp.mpf('0.5')
    root=mp.findroot(rightY_mp, (guess,guess+mp.mpf('0.1')))
    lo,hi=bracket_root(root)
    kI=sign_change(rightY_iv,lo,hi,f"R{j}")
    assert kI.a > (2*j-1)*iv.pi and kI.b < 2*j*iv.pi
    right_ks[j]=kI

YLc=leftY_iv(c_k); YRc=rightY_iv(c_k)
C0=(YRc/(1-c/R1))/(YLc/(1-c/L1))
assert C0.a>0
A={}; rho={}
for j in range(2,N+1):
    k=right_ks[j]; Rj=k*k
    Aj=-2*C0*k*leftY_iv(k)/rightYp_iv(k)*(1-Rj/R1)/(1-Rj/L1)
    A[j]=Aj
    rho[j]=Rj-c
    assert rho[j].a>0
    assert rho[j].a > (iv.pi**2*j*j).b

k2=right_ks[2]
ML=11+11*kap/(5*k2)
B=8*iv.pi*C0*ML*iv.sqrt(1+(kap/k2)**2)
for j in range(2,N+1):
    assert max(abs(A[j].a),abs(A[j].b)) < (B*j).a

pairs=[]
for i in range(2,N+1):
    for j in range(i+1,N+1):
        lam=rho[i]+rho[j]
        coeff=A[i]*A[j]*(rho[j]-rho[i])**2
        pairs.append((coeff,lam))

def q_and_qp(tI):
    q=I(0); qp=I(0)
    for coeff,lam in pairs:
        term=coeff*iv.exp(-lam*tI)
        q += term
        qp -= lam*term
    return q,qp

def U_bound(m,tau,Ntail=None):
    begin=2 if Ntail is None else Ntail+1
    JMAX=max(200,begin+200)
    total=I(0)
    for j in range(begin,JMAX+1):
        total += j**(2*m+1)*iv.exp(-(iv.pi**2)*j*j*tau)
    j=JMAX+1
    term=j**(2*m+1)*iv.exp(-(iv.pi**2)*j*j*tau)
    ratio=((j+1)/j)**(2*m+1)*iv.exp(-(iv.pi**2)*(2*j+1)*tau)
    assert ratio.b < 1
    total += term/(1-ratio)
    return B*(4*iv.pi**2)**m*total

def tail_bounds(tau,Ntail):
    U0=U_bound(0,tau); U2=U_bound(2,tau)
    U0N=U_bound(0,tau,Ntail); U2N=U_bound(2,tau,Ntail)
    return 2*(U0*U2N+U2*U0N)

# Strict negative sign on [0.0015,0.010].
tau=I('0.0015')
EQ=tail_bounds(tau,N)
def second_derivative_abs_bound(tI):
    total=I(0)
    for coeff,lam in pairs:
        total += abs(coeff)*lam**2*iv.exp(-lam*tI)
    return total

MBOX=300
lo0=mp.mpf('0.0015'); hi0=mp.mpf('0.010')
worst_upper=None
for n in range(MBOX):
    lo=lo0+(hi0-lo0)*n/MBOX
    hi=lo0+(hi0-lo0)*(n+1)/MBOX
    loI=I(mp.nstr(lo,45))
    q,qp=q_and_qp(loI)
    M2=second_derivative_abs_bound(loI)
    width=hi-lo
    end_upper=q.b+qp.b*width+I('0.5')*M2.b*width**2
    upper=max(q.b,end_upper)
    if worst_upper is None or upper>worst_upper:
        worst_upper=upper
assert worst_upper+EQ.b < 0

# First-pair dominance for all t>=0.012.
C23=A[2]*A[3]*(rho[3]-rho[2])**2
lam23=rho[2]+rho[3]
assert C23.a > 0
for j in range(3,N+1):
    assert rho[j].a > rho[j-1].b
late_tau=I('0.012')
EQ30=tail_bounds(late_tau,30)
Efinite=I(0)
for i in range(2,31):
    for j in range(i+1,31):
        if (i,j)==(2,3):
            continue
        coeff=A[i]*A[j]*(rho[j]-rho[i])**2
        lam=rho[i]+rho[j]
        Efinite += abs(coeff)*iv.exp(-(lam-lam23)*late_tau)/C23
Etail=EQ30*iv.exp(lam23*late_tau)/C23
Etotal=Efinite+Etail
assert Etotal.b < 1

print("ROBIN_Q50_NEGATIVE_AND_LATE_POSITIVE: PASS")
print("geometry: w1=25/66,D1=25; w2=4/33,D2=1/4; right width=1/2,D=1;kappa=4")
print("cutoff_N",N,"Taylor_boxes",MBOX,"interval_digits",iv.dps)
print("C0",C0)
print("c",c)
print("B",B)
print("Q_tail_error_at_tau_0.0015",EQ)
print("Q50_global_Taylor_upper",worst_upper)
print("QZ_upper_on_0.0015_to_0.010",worst_upper+EQ.b)
print("C23",C23)
print("lambda23",lam23)
print("Q_tail_error_N30_at_tau_0.012",EQ30)
print("E_finite_pairs_N30",Efinite)
print("E_tail_pairs_gt30",Etail)
print("E_total",Etotal)
print("runtime_seconds",time.time()-start)
print("environment",platform.platform(),sys.version.split()[0],mp.__version__,"machine",platform.machine(),"processor",platform.processor(),"cpus",os.cpu_count())
