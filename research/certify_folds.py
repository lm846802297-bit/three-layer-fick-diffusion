"""Interval certificates for two time-double-zeros and the steady-tie tail.
Run: OPENBLAS_NUM_THREADS=1 python certify_folds.py
"""
from fold import make
import mpmath as mp
import numpy as np
M,I=mp.mp,mp.iv
M.dps=67; I.dps=61
N=32
def s(x): return M.nstr(x,65)
def calc(k,K,z,ctx,full=False):
 c=ctx.mpf; y=c(1);p=K/k; norm=c(0); nodes=[]
 for i,(h,u) in enumerate(((c('.5')-z,c('.04')),(z,c('1.1')),(c('.5'),c(1)))):
  a,b,x=y,u*p,h*u*k
  co,si=ctx.cos(x),ctx.sin(x)
  if full:
   norm+=(a*a+b*b)*h/2+(a*a-b*b)*ctx.sin(2*x)/(4*u*k)+a*b*si**2/(u*k)
   arg=ctx.atan2(b,a)
   nodes.append(((x-arg-ctx.pi/2)/ctx.pi,(-arg-ctx.pi/2)/ctx.pi))
  y,p=a*co+b*si,(b*co-a*si)/u
  if i==1:mid=y
 if full:return K*mid*(1-y)/norm/k**2,nodes
 return p+K/k*y
def z(ctx,tie):
 c=ctx.mpf
 return c('.5')*(1-c('.04')**2)/(c('1.1')**2-c('.04')**2) if tie else c('.4132')
def roots(aa,bb,tie=False,n=N):
 a,b=map(M.mpf,(aa,bb));m=(a+b)/2
 K=I.mpf([aa,bb]); Z=z(I,tie)
 guesses=make(float(m),z=float(z(M,tie)),n=n)[1]
 L=[]; C=[]; prev=0
 def ff(x):
  l,r=M.mpf(x.a),M.mpf(x.b)
  assert M.floor(l)==M.floor(r) and l>M.floor(l)
  return int(M.floor(l))
 for g in np.sqrt(guesses):
  def single(p):
   r=M.findroot(lambda k:calc(k,p,z(M,tie),M),(M.mpf(float(g))-.001,M.mpf(float(g))+.001))
   d=M.mpf('1e-35')
   lo,hi=s(r-d),s(r+d)
   assert calc(I.mpf(lo),I.mpf(s(p)),Z,I)*calc(I.mpf(hi),I.mpf(s(p)),Z,I)<0
   return lo,hi
  lo=single(a)[0]; hi=single(b)[1]
  # Ordered mode's eigenvalue is nondecreasing in kappa by Robin minmax.
  # The interval node counts below identify both endpoint roots and the whole branch.
  k=I.mpf([lo,hi]);assert k.a>prev;prev=k.b
  w,ns=calc(k,K,Z,I,True)
  assert sum(ff(t)-ff(v) for t,v in ns)==len(L)
  L.append(k**2);C.append(w)
 return L,C
def val(t,d,L,C,K):
 # |flux coefficient| <= 17.34*(1+2.2*sqrt(lambda)).
 # Past mode 32, lambda>=8*j^2. For j>=32, t>=.023,
 # exp(-lambda*t)*lambda^p decreases in lambda for p<=1.5.
 # At lambda=8*j^2 the bound for all orders <= 2 is <=1000*j^3.
 # Its ratio is <=(33/32)^3*exp(-8*65*.023)<1/2.
 assert t.a>I.mpf('.023') and len(L)==N
 tail=2000*N**3*I.exp(-8*N**2*t.a)
 rr=I.mpf('.5')+1/K; Z=z(I,False)
 rl=Z*I.mpf('1.1')**2+(I.mpf('.5')-Z)*I.mpf('.04')**2+1/K
 return ((rr-rl)/(rr+rl) if d==0 else 0)-sum(c*(-l)**d*I.exp(-l*t) for l,c in zip(L,C))+I.mpf([-1,1])*tail
def fold(ka,kb,ta,tb,sg):
 ta,tb=map(M.mpf,(ta,tb));tc=(ta+tb)/2
 a,b=map(M.mpf,(ka,kb))
 def tester(k1,k2):
  k1,k2=s(k1),s(k2); L,C=roots(k1,k2); K=I.mpf([k1,k2])
  def ev(p,q,d):return val(I.mpf([s(p),s(q)]),d,L,C,K)
  return ev
 for l,r in zip(np.linspace(0,1,65)[:-1],np.linspace(0,1,65)[1:]):
  ev=tester(a+M.mpf(float(l))*(b-a),a+M.mpf(float(r))*(b-a))
  assert sg*ev(ta,ta,1)<0 and sg*ev(tb,tb,1)>0,(ev(ta,ta,1),ev(tb,tb,1))
  assert sg*ev(ta,tb,2)>0,ev(ta,tb,2)
 for k,j in [(a,-1),(b,1)]:
  ev=tester(k,k)
  vc=ev(tc,tc,0); dy=ev(tc,tc,1); ddy=ev(ta,tb,2)
  r=I.mpf(s((tb-ta)/2));B=abs(dy)*r+abs(ddy)*r*r/2
  assert j*vc-B>0,(vc,B)
  print('PASS sign on entire time window',s(k),j,'central value',vc,'Taylor margin',B,flush=True)
 print('PASS strict curvature and time-edge derivative signs:',ka,kb,s(ta),s(tb),sg,flush=True)
if __name__=='__main__':
 fold('5.79250','5.79255','.25470','.25481',-1)
 fold('8.66518','8.66522','.023891','.023918',1)
 for k,sg in [('1',1),('10',-1)]:
  L,C=roots(k,k,tie=True,n=1)
  assert sg*C[0]>0
  print('PASS tie: kappa=',k,'lambda_1=',L[0],'c_1=',C[0],flush=True)
 # A principal-mode zero on steady tie; the second mode remains positive
 # throughout this interval, so at that zero the late sign is negative.
 ka,kb='5.69731053','5.69731055'
 _,ca=roots(ka,ka,tie=True,n=1)
 _,cb=roots(kb,kb,tie=True,n=1)
 L,C=roots(ka,kb,tie=True,n=2)
 assert ca[0]>0 and cb[0]<0 and C[1]>0
 print('PASS tie principal-mode zero:',ka,kb,'lambda_2=',L[1],'c_2=',C[1],flush=True)
