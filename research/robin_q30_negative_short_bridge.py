import platform, sys
import mpmath as mp
from mpmath import iv
mp.mp.dps=70
iv.dps=60
I=iv.mpf
kap=I(4); w1=I(25)/66; w2=I(4)/33; wr=I(1)/2
D1=I(25); D2=I(1)/4; DR=I(1)
def mat_iv(k,h,D):
 z=k*h/iv.sqrt(D)
 return [[iv.cos(z),iv.sin(z)/(iv.sqrt(D)*k)],[-iv.sqrt(D)*k*iv.sin(z),iv.cos(z)]]
def mm(A,B):
 return [[A[0][0]*B[0][0]+A[0][1]*B[1][0],A[0][0]*B[0][1]+A[0][1]*B[1][1]],
         [A[1][0]*B[0][0]+A[1][1]*B[1][0],A[1][0]*B[0][1]+A[1][1]*B[1][1]]]
def leftY_iv(k):
 M=mm(mat_iv(k,w2,D2),mat_iv(k,w1,D1)); return M[0][0]+M[0][1]*kap
def left_phase_iv(k):
 theta=k*I(5)/66-iv.atan2(I(4)/5,k)
 return iv.atan2(10*iv.sin(theta),iv.cos(theta))+2*k*w2
def rightY_iv(k): return iv.cos(k/2)+kap/k*iv.sin(k/2)
def rightYp_iv(k): return -I(1)/2*iv.sin(k/2)+kap*(I(1)/2*iv.cos(k/2)/k-iv.sin(k/2)/k**2)
def fullF_iv(k):
 M=mm(mat_iv(k,wr,DR),mm(mat_iv(k,w2,D2),mat_iv(k,w1,D1)))
 y=M[0][0]+M[0][1]*kap; p=M[1][0]+M[1][1]*kap
 return p+kap*y
def mat_mp(k,h,D):
 z=k*h/mp.sqrt(D); c=mp.cos(z); s=mp.sin(z)
 return [[c,s/(mp.sqrt(D)*k)],[-mp.sqrt(D)*k*s,c]]
def mm_mp(A,B):
 return [[sum(A[i][q]*B[q][j] for q in range(2)) for j in range(2)] for i in range(2)]
def leftY_mp(k):
 M=mm_mp(mat_mp(k,mp.mpf(4)/33,mp.mpf(1)/4),mat_mp(k,mp.mpf(25)/66,25)); return M[0][0]+M[0][1]*4
def rightY_mp(k): return mp.cos(k/2)+4/k*mp.sin(k/2)
def fullF_mp(k):
 M=mm_mp(mat_mp(k,mp.mpf(1)/2,1),mm_mp(mat_mp(k,mp.mpf(4)/33,mp.mpf(1)/4),mat_mp(k,mp.mpf(25)/66,25)))
 y=M[0][0]+4*M[0][1]; p=M[1][0]+4*M[1][1]
 return p+4*y
def bracket(root):
 rad=mp.mpf('1e-25')
 return mp.nstr(root-rad,65),mp.nstr(root+rad,65)
def root_iv(fun,root,label):
 lo,hi=bracket(root); fl=fun(I(lo)); fh=fun(I(hi))
 assert (fl.b<0 and fh.a>0) or (fl.a>0 and fh.b<0),(label,fl,fh)
 return I([lo,hi])
kL1mp=mp.findroot(leftY_mp,(3.7,3.8)); kR1mp=mp.findroot(rightY_mp,(4.5,4.7)); kcmp=mp.findroot(fullF_mp,(3.9,4.1))
kL1=root_iv(leftY_iv,kL1mp,'L1'); kR1=root_iv(rightY_iv,kR1mp,'R1'); kc=root_iv(fullF_iv,kcmp,'lambda2')
phase0=left_phase_iv(kL1.a); phase1=left_phase_iv(kL1.b)
assert phase0.b<iv.pi/2 and phase1.a>iv.pi/2
th0=kL1.a*I(5)/66-iv.atan2(I(4)/5,kL1.a); th1=kL1.b*I(5)/66-iv.atan2(I(4)/5,kL1.b)
assert th0.a>-iv.pi/2 and th1.b<iv.pi/2
L1=kL1*kL1; R1=kR1*kR1; c=kc*kc
assert L1.a<c.a and c.b<R1.a
right={1:kR1}
for j in range(2,71):
 guess=(2*j-1)*mp.pi+mp.mpf('0.5')
 kr=mp.findroot(rightY_mp,(guess,guess+mp.mpf('0.1')))
 kj=root_iv(rightY_iv,kr,'R'+str(j))
 assert kj.a>(2*j-1)*iv.pi and kj.b<2*j*iv.pi
 right[j]=kj
YLc=leftY_iv(kc); YRc=rightY_iv(kc)
C0=(YRc/(1-c/R1))/(YLc/(1-c/L1))
A={}; rho={}
for j in range(2,71):
 k=right[j]; Rj=k*k
 A[j]=-2*C0*k*leftY_iv(k)/rightYp_iv(k)*(1-Rj/R1)/(1-Rj/L1)
 rho[j]=Rj-c
 assert rho[j].a>0
pairs30=[]
for i in range(2,31):
 for j in range(i+1,31):
  pairs30.append((A[i]*A[j]*(rho[j]-rho[i])**2,rho[i]+rho[j]))
def q30(t):
 q=I(0)
 for coeff,lam in pairs30:
  q+=coeff*iv.exp(-lam*t)
 return q
qleft=q30(I('0.0015')); qright=q30(I('0.004'))
M=2000
max_upper=None; min_lower=None
for n in range(M):
 lo=mp.mpf('0.0015')+mp.mpf('0.0025')*n/M
 hi=mp.mpf('0.0015')+mp.mpf('0.0025')*(n+1)/M
 qi=q30(I([mp.nstr(lo,45),mp.nstr(hi,45)]))
 max_upper=qi.b if max_upper is None else max(max_upper,qi.b)
 min_lower=qi.a if min_lower is None else min(min_lower,qi.a)
T30_70=I(0)
for i in range(2,71):
 for j in range(i+1,71):
  if j>30:
   cij=A[i]*A[j]*(rho[j]-rho[i])**2
   T30_70+=abs(cij)*iv.exp(-(rho[i]+rho[j])*I('0.001'))
tau=I('0.001'); k2=right[2]
ML=11+11*kap/(5*k2)
B=8*iv.pi*C0*ML*iv.sqrt(1+(kap/k2)**2)
def U(m,N=None):
 start=2 if N is None else N+1; J=240; total=I(0)
 for j in range(start,J+1): total+=j**(2*m+1)*iv.exp(-iv.pi**2*j*j*tau)
 j=J+1; term=j**(2*m+1)*iv.exp(-iv.pi**2*j*j*tau)
 ratio=((j+1)/j)**(2*m+1)*iv.exp(-iv.pi**2*(2*j+1)*tau)
 assert ratio.b<1
 return B*(4*iv.pi**2)**m*(total+term/(1-ratio))
U0=U(0); U2=U(2); UN0=U(0,70); UN2=U(2,70)
Qtail70=2*(U0*UN2+U2*UN0)
Qtail30=T30_70+Qtail70
QZupper=max_upper+Qtail30.b
assert max_upper<0 and QZupper<0
print('ROBIN_Q30_NEGATIVE_SHORT_BRIDGE: PASS')
print('cutoff N=30; root certificates through R70; time window [0.0015,0.004]; interval boxes=2000')
print('Q30(0.0015)',qleft); print('Q30(0.004)',qright)
print('Q30 maximum upper bound',max_upper)
print('Q30 minimum lower bound',min_lower)
print('B',B)
print('finite omitted pairs 30<j<=70 at tau=0.001',T30_70)
print('Gaussian omitted tail j>70 at tau=0.001',Qtail70)
print('total sup |Q-Q30| for t>=0.001',Qtail30)
print('sup QZ on [0.0015,0.004]',QZupper)
print('runtime',platform.platform(),'Python',sys.version.split()[0],'mpmath',mp.__version__)
