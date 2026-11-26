"""Certified spectral slabs for the fixed fast-at-wall example."""
import mpmath as mp
import numpy as np
from spectral import spectrum
mp.mp.dps=76; mp.iv.dps=67
M=mp.mp; I=mp.iv
mpf=M.mpf; iv=I.mpf
K='6.5'
layers=[('.0868','.04'),('.4132','1.1'),('.5','1')]
def calc(k,ctx,full=False):
 c=ctx.mpf; kap=c(K);y=c(1);p=kap/k
 norm=c(0); nodes=[]
 for i,(hh,uu) in enumerate(layers):
  h,u=c(hh),c(uu)
  a,b,x=y,u*p,h*u*k
  co,si=ctx.cos(x),ctx.sin(x)
  if full:
   norm+=(a*a+b*b)*h/2+(a*a-b*b)*ctx.sin(2*x)/(4*u*k)+a*b*(si**2)/(u*k)
   arg=ctx.atan2(b,a)
   nodes.append(((x-arg-ctx.pi/2)/ctx.pi,(-arg-ctx.pi/2)/ctx.pi))
  y,p=a*co+b*si,(b*co-a*si)/u
  if i==1:mid=y
 if full:return p+kap/k*y,kap*mid*(1-y)/norm/k**2,nodes
 return p+kap/k*y

def cert_roots(N):
 # These floats only locate candidate brackets; every bracket and the count are certified below.
 guesses,_=spectrum(N)
 L=[]; C=[]; prev=0
 def ff(x):
  l=mpf(x.a); r=mpf(x.b)
  assert M.floor(l)==M.floor(r) and l>M.floor(l)
  return int(M.floor(l))
 for g in np.sqrt(guesses):
  r=M.findroot(lambda k:calc(k,M),(mpf(float(g))-.001,mpf(float(g))+.001),solver='secant')
  lo,hi=[M.nstr(r+d*mpf('1e-48'),73) for d in (-1,1)]
  a,b=iv(lo),iv(hi)
  f,g=calc(a,I),calc(b,I)
  assert f*g<0 and a>prev
  prev=b
  k=iv([lo,hi]); _,weight,nodes=calc(k,I,True)
  n=sum(ff(t)-ff(s) for t,s in nodes)
  assert n==len(L),(n,len(L))
  L.append(k*k); C.append(weight)
 return L,C

# Neumann comparison lambda_n >= Dmin*pi^2*(n-1)^2 > 8*(n-1)^2,
# n starts at 1. For normalized eigenfunctions, |phi|^2 <= 1+2*sqrt(lambda/Dmin).
# Therefore |flux coeff| <= 13*(1+2.2*sqrt(lambda)).
# For n>N and t>=.000999, 1 and sqrt(lambda) times exp(-lambda*t)
# are decreasing in lambda. a_n <= 100*n*exp(-8*n*n*t), n=N,N+1,...
# Ratio <= 2 exp(-8*(2N+1)*.000999) < 1/2 for N=120.
# Same bound also works for integrated coeff because omitted lambda > 1.
def tail(N,left):
 return 200*N*I.exp(-8*N*N*left)

def certify_early():
 from fractions import Fraction as F
 from math import factorial
 def E(x):return sum((x**j/F(factorial(j)) for j in range(31)),F(0))
 d1=F('.0868')*F('.04'); d2=F('.4132')*F('1.1'); dR=F('.5')
 dL=d1+d2; kap=F('6.5'); T=F('.001'); z=dL/(2*T)
 qR=dR/(2*T); gamma=(1+25/F(10,11))*(1+F(10,11))/8
 C=kap/(2*gamma); g=(dR**2-dL**2)/4
 U=[1/E(2*d*z) for d in (d1,d2,dR)]
 EW=(1+U[0])*(1+U[1])*(1+U[2])-1
 assert EW<F('.22') and U[2]<F('.001')
 # sqrt(pi)>1.73, sqrt(T)<.032.
 assert (EW+U[2])/(1-EW)+4*F('.032')/(dL*F('1.73'))+kap*F('.04')/z<F('.5')
 assert 1/E(2*dR*qR)<F(1,3)
 assert g>T/2
 # sqrt(pi)<1.8, sqrt(T)>.0316. The ratio bound is increasing up to T.
 assert 3*kap*dL**2*F('1.8')/(C*dR*F('.0316')*E(g/T))<F('.06')
 print('PASS 0<t<=.001, right/left < .06 using contour and Chernoff bounds')

def certify(L,C):
 N=len(L); f0=-iv(9009)/106259009
 # Each item: (time slab, order (0=value, 1=derivative), sign).
 plan=[('.000999','.012',0,1),('.012','.015',1,-1),('.015','.08',0,-1),
       ('.08','.09',1,1),('.09','.55',0,1),('.55','.65',1,-1),('.65','1',0,-1)]
 # Mutable stats reassure that recursion verified all slabs.
 for aa,bb,order,sg in plan:
  a,b=mpf(aa),mpf(bb)
  def val(t,order=0):
   r=sum((c*((-l)**order)*I.exp(-l*iv(t)) for l,c in zip(L,C)),iv(0))
   return (f0 if order==0 else iv(0))-r
  if order:
   assert val(a)*sg+tail(N,iv(a))<0 and val(b)*sg-tail(N,iv(b))>0,(aa,bb,val(a),val(b),sg)
  tests=0; maxdepth=0; accepted=0
  P=12
  weights=[(-c)*((-l)**order) for l,c in zip(L,C)]
  def test(a,b,d):
   nonlocal tests,maxdepth,accepted
   tests+=1;maxdepth=max(maxdepth,d)
   t=iv([a,b]); mid=(iv(a)+iv(b))/2; rad=(iv(b)-iv(a))/2
   Q=[iv(0) for _ in range(P)]
   if order==0:Q[0]=f0
   rem=iv(0)
   for l,w in zip(L,weights):
    f=w*I.exp(-l*mid)
    for j in range(P):
     Q[j]+=f
     f=f*(-l)/(j+1)
    g=abs(w)*I.exp(-l*iv(a))
    rem+=g*(l*rad)**P/M.factorial(P)
   span=sum((abs(Q[j])*rad**j for j in range(1,P)),iv(0))+rem+tail(N,iv(a))
   if Q[0]*sg-span>0:
    accepted+=1;return
   assert d<28,(aa,bb,a,b,Q[0],span,order)
   t=(a+b)/2
   test(a,t,d+1);test(t,b,d+1)
  test(a,b,0)
  print('PASS',aa,bb,'derivative' if order else 'value',sg,'cells',accepted,'depth',maxdepth,flush=True)
 for a,b,sg in [('.013569414','.013569416',-1),('.083737507','.083737510',1),('.596764946','.596764949',-1)]:
  a,b=mpf(a),mpf(b)
  assert val(a)*sg+tail(N,iv(a))<0 and val(b)*sg-tail(N,iv(b))>0
  print('PASS refined root bracket',mp.nstr(a,12),mp.nstr(b,12))
 late=sum((abs(c)*I.exp(-l) for l,c in zip(L,C)),iv(0))+tail(N,iv(1))
 assert late<-f0
 print('PASS t>=1',str(late),str(-f0))

if __name__=='__main__':
 L,C=cert_roots(120);print('PASS 120 isolated modes, ordered by node count',flush=True)
 certify_early()
 certify(L,C)
