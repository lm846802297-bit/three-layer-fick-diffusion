import json, math, platform, sys, subprocess
import mpmath as mp
from mpmath import iv

mp.mp.dps = 90
iv.dps = 70
I = iv.mpf
kap = I(4)
w1, w2, wr = I(25)/66, I(4)/33, I(1)/2
D1, D2, DR = I(25), I(1)/4, I(1)

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

# High precision scalar versions for proposing narrow brackets.
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

# First left root, first right root, and global lambda_2.
kL1_mp=mp.findroot(leftY_mp, (3.7,3.8))
kR1_mp=mp.findroot(rightY_mp, (4.5,4.7))
kG2_mp=mp.findroot(fullF_mp, (3.9,4.1))
L1_lo,L1_hi=bracket_root(kL1_mp)
L1k=sign_change(leftY_iv,L1_lo,L1_hi,"L1")
# The left Prüfer phase starts at zero, is strictly increasing, and reaches pi/2
# at the first left Dirichlet eigenvalue; this certifies the first-mode index.
phase_lo=left_phase_iv(I(L1_lo)); phase_hi=left_phase_iv(I(L1_hi))
assert phase_lo.b < iv.pi/2 and phase_hi.a > iv.pi/2
theta_lo=I(L1_lo)*I(5)/66-iv.atan2(I(4)/5,I(L1_lo))
theta_hi=I(L1_hi)*I(5)/66-iv.atan2(I(4)/5,I(L1_hi))
assert theta_lo.a > -iv.pi/2 and theta_hi.b < iv.pi/2
R1k=sign_change(rightY_iv,*bracket_root(kR1_mp),"R1")
c_k=sign_change(fullF_iv,*bracket_root(kG2_mp),"global lambda2")
L1=L1k**2; R1=R1k**2; c=c_k**2
assert L1.a < c.a and c.b < R1.a
# Midpoint Dirichlet interlacing identifies the unique global eigenvalue in (L1,R1) as lambda_2.

# Right roots R_j, j=2,...,30. Their Prüfer phase is strictly increasing:
# theta_R'(k)=1/2+4/(k^2+16)>0, and root j lies in ((2j-1)pi,2j*pi).
right_ks={1:R1k}
for j in range(2,31):
    guess=(2*j-1)*mp.pi + mp.mpf('0.5')
    root=mp.findroot(rightY_mp, (guess,guess+mp.mpf('0.1')))
    lo,hi=bracket_root(root)
    kI=sign_change(rightY_iv,lo,hi,f"R{j}")
    assert kI.a > (2*j-1)*iv.pi and kI.b < 2*j*iv.pi
    right_ks[j]=kI

# C0 normalizes the tail quotient so Phi_Z(0)=1.
YLc=leftY_iv(c_k); YRc=rightY_iv(c_k)
C0=(YRc/(1-c/R1))/(YLc/(1-c/L1))
assert C0.a>0

# Residues A_j for j=2,...,30. This formula removes the distinct first roots L1 and R1.
A={}; rho={}
for j in range(2,31):
    k=right_ks[j]; Rj=k*k
    Aj=-2*C0*k*leftY_iv(k)/rightYp_iv(k)*(1-Rj/R1)/(1-Rj/L1)
    A[j]=Aj
    rho[j]=Rj-c
    assert rho[j].a>0

# Check |A_j| <= B*j, using the analytic right-root derivative bound and the
# fact 0 < (1-R_j/R_1)/(1-R_j/L_1) < 1 for j>=2.
k2=right_ks[2]
ML=11+11*kap/(5*k2)
B=8*iv.pi*C0*ML*iv.sqrt(1+(kap/k2)**2)
for j in range(2,31):
    assert max(abs(A[j].a),abs(A[j].b)) < (B*j).a

# Build Q_30 and Q_30' as interval sums.
pairs=[]
for i in range(2,31):
    for j in range(i+1,31):
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

T0=I('0.010'); T1=I('0.012')
q0,_=q_and_qp(T0); q1,_=q_and_qp(T1)
qroot_lo,_=q_and_qp(I('0.01102'))
qroot_hi,_=q_and_qp(I('0.01103'))
assert q0.b < 0 and q1.a > 0
assert qroot_lo.b < 0 and qroot_hi.a > 0

# Certify Q_30'>0 on [0.010,0.012] with 2000 interval sub-boxes.
M=2000
qp_lowers=[]; qp_uppers=[]
for n in range(M):
    lo=mp.mpf('0.010')+mp.mpf('0.002')*n/M
    hi=mp.mpf('0.010')+mp.mpf('0.002')*(n+1)/M
    _,qpi=q_and_qp(I([mp.nstr(lo,40),mp.nstr(hi,40)]))
    qp_lowers.append(qpi.a); qp_uppers.append(qpi.b)
qp_lower=min(qp_lowers); qp_upper=max(qp_uppers)
assert qp_lower>0

# Gaussian tail bounds. For each m, sum to JMAX then use a geometric tail.
def U_bound(m,N=None):
    start=2 if N is None else N+1
    JMAX=160
    total=I(0)
    for j in range(start,JMAX+1):
        total += j**(2*m+1)*iv.exp(-(iv.pi**2)*j*j*T0)
    j=JMAX+1
    term=j**(2*m+1)*iv.exp(-(iv.pi**2)*j*j*T0)
    ratio=((j+1)/j)**(2*m+1)*iv.exp(-(iv.pi**2)*(2*j+1)*T0)
    assert ratio.b < 1
    total += term/(1-ratio)
    return B*(4*iv.pi**2)**m*total
U=[U_bound(m) for m in range(4)]
UN=[U_bound(m,30) for m in range(4)]
EQ=2*(U[0]*UN[2]+U[2]*UN[0])
EQP=2*(U[3]*UN[0]+U[2]*UN[1]+U[1]*UN[2]+U[0]*UN[3])
assert qp_lower > EQP.b
assert q0.b + EQ.b < 0
assert q1.a - EQ.b > 0
assert qroot_lo.b + EQ.b < 0
assert qroot_hi.a - EQ.b > 0

print("ROBIN_Q30_MIDDLE_WINDOW_CERTIFICATE: PASS")
print("geometry: w1=25/66,D1=25; w2=4/33,D2=1/4; right length=1/2,D=1;kappa=4")
print("root brackets:")
print("  k_L1",L1k)
print("  k_R1",R1k)
print("  k_lambda2",c_k)
print("  k_R2",right_ks[2])
print("  k_R30",right_ks[30])
print("c=lambda2",c)
print("C0",C0)
print("B",B)
print("Q30(0.010)",q0)
print("Q30(0.012)",q1)
print("Q30(0.01102)",qroot_lo)
print("Q30(0.01103)",qroot_hi)
print("Q30_prime lower bound over 2000 subintervals",qp_lower)
print("Q30_prime upper bound over 2000 subintervals",qp_upper)
print("U_m",U)
print("U_m,N30",UN)
print("sup |Q-Q30|",EQ)
print("sup |Q'-Q30'|",EQP)
print("Q_Z derivative lower bound",qp_lower-EQP.b)
print("Q_Z(0.010) upper bound",q0.b+EQ.b)
print("Q_Z(0.012) lower bound",q1.a-EQ.b)
print("Q_Z(0.01102) upper bound",qroot_lo.b+EQ.b)
print("Q_Z(0.01103) lower bound",qroot_hi.a-EQ.b)
print("runtime",platform.platform(),"Python",sys.version.split()[0],"mpmath",mp.__version__)
