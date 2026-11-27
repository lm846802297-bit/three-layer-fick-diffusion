"""Uniform 3-layer static-tie cross-check with outward mpmath intervals.
Rational Robin endpoints, outwardly checked root proposals; no sampled-time decisions.
Python >=3.9, scipy (only supplies proposals), mpmath 1.3.0. Run directly.
"""
from fractions import Fraction as F
from functools import lru_cache
import math as M
from scipy.optimize import brentq
from mpmath import iv
import sys
iv.dps=35
I=iv.mpf; P=iv.pi
l1=I(5)/66;l2=I(8)/33;d=I(2)/11;u=I(2);v=I(1)/5
N=25
def rat(q): return I(q.numerator)/q.denominator
def join(a,b):return I([a.a,b.b])
def YL(k,K):
 c1,c2=iv.cos(k*l1),iv.cos(k*l2);s1,s2=iv.sin(k*l1),iv.sin(k*l2)
 return (c1+K*v/k*s1)*c2+(K*u/k*c1-u/v*s1)*s2
def YR(k,K):return iv.cos(k/2)+K/k*iv.sin(k/2)
def phase(k,K):
 d0=k*l1-iv.atan2(K*v,k)
 # Left L1/L2 and right R1/R2/R3 lie in this branch for this geometry.
 assert -P/2<d0<P/2
 return iv.atan2(u/v*iv.sin(d0),iv.cos(d0))+k*l2
def full(k,K):
 c1,c2=iv.cos(k*l1),iv.cos(k*l2); s1,s2=iv.sin(k*l1),iv.sin(k*l2)
 yl=YL(k,K); yr=YR(k,K)
 y0=c1+K*v/k*s1;q0=K*c1-k/v*s1
 return yl*(K*iv.cos(k/2)-k*iv.sin(k/2))+yr*(q0*c2-k/u*y0*s2)
def cert_root(fun,x,inc=True): # integer grid: proposals, not evidence
 n=M.floor(x*10**10); l=rat(F(n-1,10**10));r=rat(F(n+2,10**10))
 assert (fun(l)<0<fun(r)) if inc else (fun(r)<0<fun(l))
 return join(l,r)
@lru_cache(None)
def roots(KF):
 K=rat(KF); kf=float(KF)
 def ph(k):
  D=k*5/66-M.atan2(kf*.2,k)
  j=M.floor((D+M.pi/2)/M.pi); D-=j*M.pi
  return j*M.pi+M.atan2(10*M.sin(D),M.cos(D))+k*8/33
 k1=brentq(lambda k:ph(k)-M.pi/2,1e-10,20)
 A=cert_root(lambda k:phase(k,K)-P/2,k1)
 rr=[]
 for n in range(1,N+1):
  u0=(2*n-1)*M.pi
  x=brentq(lambda k:k-u0-2*M.atan2(kf,k),u0-1e-8,2*n*M.pi+1e-8)
  rr.append(cert_root(lambda k:k-(2*n-1)*P-2*iv.atan2(K,k),x))
 if kf==0:return A,rr,I(0)
 def nf(k):
  l=k*5/66;r=k*8/33
  c1,c2,s1,s2=M.cos(l),M.cos(r),M.sin(l),M.sin(r)
  y=c1+kf*.2/k*s1; q=kf*c1-5*k*s1
  yl=y*c2+2*q/k*s2; yr=M.cos(k/2)+kf/k*M.sin(k/2)
  return yl*(kf*M.cos(k/2)-k*M.sin(k/2))+yr*(q*c2-.5*k*y*s2)
 grd=brentq(nf,1e-12,min(k1,float(rr[0].a))-1e-9)
 G=cert_root(lambda k:full(k,K),grd,False)
 assert G<A and G<rr[0]
 return A,rr,G

def tail_const(n,t,j): # B is for Htilde=(b/a)H
 p=1+2*j
 first=100*(2*P*n)**p*iv.exp(-((2*n-1)*P)**2*t)
 q=(I(n+1)/n)**p*iv.exp(-8*P*P*n*t)
 assert q<1
 return (first/(1-q)).b
def prepare(l,r):
 K=join(rat(l),rat(r));a1,u1,g1=roots(l);a2,u2,g2=roots(r)
 a=join(a1,a2)**2;ks=[join(v,w) for v,w in zip(u1,u2)]; b=ks[0]**2;m=ks[1]**2
 if not (a<b and all((phase(k,K)>(P/2) if n==1 else phase(k,K)<(n-I(1)/2)*P) for n,k in enumerate(ks[:3],1))):return
 # n>=4: continuous positive-frequency interface phase gain <pi/2,
 # and the right Robin endpoint adds <pi/2 to its pole wavenumber.
 # Left phase < k*(l1+l2)+pi/2; right angle k/2-atan(K/k).
 # Since ell gap 2/11 and k>(2n-1)*pi, k*gap>pi.
 assert 7*d>1
 # Unscaled residues satisfy |Btilde| <= C*k_n for n>=2.
 assert 4*(11+K*(u+v)/(3*P))*iv.sqrt(1+(K/(3*P))**2)<100
 # For all times up to t0, fixed bound B0 and positivity of S.
 t0=I('0.0008');Q=d/(2*t0)
 yd=((1-K/Q)-(1+K/Q)*iv.exp(-Q))/2
 mult=(1+b/Q**2)/(1-a/Q**2) # a/b<=1 on these boxes
 assert yd>0 and (11+K*(u+v)/Q)/yd*mult<30 and d**2/4>I(3)/2*t0
 B0=30/P*iv.exp(-d**2/(4*t0))*(Q*iv.sqrt(P/t0)+1/t0)
 if not B0<a*(1-t0*B0):return
 Rs=[k**2 for k in ks[1:]]
 Bs=[(1-b/R)*YL(k,K)/(iv.sin(k/2)*(1+(K**2+2*K)/R)/(4*k)) for R,k in zip(Rs,ks[1:])]
 T=I('.04')
 late=Bs[0]-sum(abs(B)*iv.exp(-(R-m)*T) for R,B in zip(Rs[1:],Bs[1:]))-iv.exp(m*T)*tail_const(N+1,T,0)
 if not late>0:return
 # Determine sign of c1 optionally. No division in the comparison.
 g=join(g1,g2)
 ratio=None if not g>0 else (YL(g,K)/YR(g,K))
 return Rs,Bs,ratio
def verify(l,r):
 D=prepare(l,r)
 if D is None:return None
 Rs,Bs,ratio=D
 def bounds(start,stop,j): # integer time unit 1e-7
  l0=I(start)/10000000;r0=I(stop)/10000000
  t=join(l0,r0)
  H=sum(B*(-R)**j*iv.exp(-R*t) for B,R in zip(Bs,Rs))
  err=tail_const(N+1,l0,j)
  return H+I([-err,err])
 def checks(left,right,ordered,fn_sign):
  pre,der,post=[],[],[]
  for L,R,j,sg in [(8000,left,0,-1),(left,right,1,1),(right,400000,0,1)]:
   q=[(L,R)]; it=0
   while q:
    ll,rr=q.pop();v=bounds(ll,rr,j)
    if sg*v>0: it+=1;continue
    # An inconclusive point cannot be fixed by subdividing the time axis.
    mm=(ll+rr)//2
    if rr-ll<=20 or not sg*bounds(mm,mm,j)>0:return None
    q.extend([(ll,mm),(mm,rr)])
   (der if j else pre if sg<0 else post).append(it)
  return pre[0],der[0],post[0]
 # a common conservative crossing bracket (or smaller ones if needed).
 for lo,hi in [(12000,18500),(11000,19000),(12500,17500),(13000,17500)]:
  ok=checks(lo,hi,0,1)
  if ok is not None:return (ok,ratio,lo,hi)
 return None
def scan(L,R,step):
 todo=[(F(i,100),F(min(i+step,R),100)) for i in reversed(range(L,R,step))];g=[];un=[]
 while todo:
  l,r=todo.pop();out=verify(l,r)
  if out:
   counts,ratio,lo,hi=out; sign='unresolved' if ratio is None else 'one' if ratio<1 else 'zero' if ratio>1 else 'unresolved'
   print(str(l),str(r),lo,hi,*counts,sign, str(ratio),flush=True);g.append((l,r,sign))
  elif r-l>F(1,500):
   m=(l+r)/2;todo += [(m,r),(l,m)]
  else:
   un.append((l,r));print('FAILED',l,r,flush=True)
 print('SUMMARY boxes',len(g),'unproven',un,flush=True)
if __name__=='__main__':
 L,R,step=map(int,sys.argv[1:] if len(sys.argv)>1 else [100,1000,50])
 scan(L,R,step)
