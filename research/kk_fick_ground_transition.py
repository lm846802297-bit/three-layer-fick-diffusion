"""Certified simple physical switch for the fixed 3-layer geometry. Run next to
kk_fick_kappa_certificate.py. Rigorous implicit differentiation using interval
automatic differentiation. No numerical differentiation or parameter sampling.
"""
from kk_fick_kappa_certificate import roots,rat,join,l1,l2,I,iv,F
from scipy.optimize import brentq
import math as M
class Jet:
 def __init__(self,x,k=0,K=0):self.x=x;self.k=k;self.K=K
 def val(o):return o if isinstance(o,Jet) else Jet(o)
 def __add__(s,o):
  o=Jet.val(o);return Jet(s.x+o.x,s.k+o.k,s.K+o.K)
 __radd__=__add__
 def __neg__(s):return Jet(-s.x,-s.k,-s.K)
 def __sub__(s,o):return s+(-Jet.val(o))
 def __rsub__(s,o):return Jet.val(o)+(-s)
 def __mul__(s,o):
  o=Jet.val(o);return Jet(s.x*o.x,s.k*o.x+s.x*o.k,s.K*o.x+s.x*o.K)
 __rmul__=__mul__
 def __truediv__(s,o):
  o=Jet.val(o);return Jet(s.x/o.x,(s.k*o.x-s.x*o.k)/o.x**2,(s.K*o.x-s.x*o.K)/o.x**2)
 def __rtruediv__(s,o):return Jet.val(o)/s
 def cos(s):
  d=-iv.sin(s.x);return Jet(iv.cos(s.x),d*s.k,d*s.K)
 def sin(s):
  d=iv.cos(s.x);return Jet(iv.sin(s.x),d*s.k,d*s.K)
def calc(lo,hi):
 _,_,ga=roots(lo);_,_,gb=roots(hi)
 k=Jet(join(ga,gb),1,0);K=Jet(join(rat(lo),rat(hi)),0,1)
 c1,c2=(k*l1).cos(),(k*l2).cos(); s1,s2=(k*l1).sin(),(k*l2).sin()
 y=c1+K/k/5*s1;ql=K*c1-k*5*s1
 yl=y*c2+ql*2/k*s2; Ql=ql*c2-k*y/2*s2
 yr=(k/2).cos()+K/k*(k/2).sin();Qr=K*(k/2).cos()-k*(k/2).sin()
 M=yl*Qr+yr*Ql; E=yl-yr
 assert yr.x>0 and M.k<0
 gE=E.K-E.k*M.K/M.k
 return E.x,gE,yl.x/yr.x
def scan(L,R,rs=10):
 chunks=[(F(i,100),F(min(i+rs,R),100)) for i in range(L,R,rs)]
 pos=[]; neg=[]; ds=[]
 while chunks:
  a,b=chunks.pop();E,g,r=calc(a,b)
  if E>0:pos.append((a,b));continue
  if E<0:neg.append((a,b));continue
  if F(6)<=a<b<=F(9) and g>0:ds.append((a,b));continue
  assert b-a>F(1,10**9),(a,b,E,g)
  mid=(a+b)/2; chunks.extend([(a,mid),(mid,b)])
 def merge(L):
  out=[]
  for a,b in sorted(L):
   if out and out[-1][1]==a:out[-1]=(out[-1][0],b)
   else:out.append((a,b))
  return [(str(a),str(b)) for a,b in out]
 n,p,der=merge(neg),merge(pos),merge(ds)
 print('E negative',n,'E positive',p,'dE/dkappa positive where sign uncertain',der)
 assert len(n)==len(p)==len(der)==1 and n[0][1]==der[0][0] and der[0][1]==p[0][0]
 
 # Exact decimal split around the transition; iterations of endpoint boxes,
 # not interval midpoint signs unless strict. No unique claim outside [1,10].
 a=F(der[0][0]);b=F(der[0][1])
 # Multiples of 1e-7 ensure clean, outward certified printed bracket.
 Lo=M.floor(a*10**7);Hi=M.ceil(b*10**7)
 while Hi-Lo>1:
  mid=(Hi+Lo)//2;E,_,_=calc(F(mid,10**7),F(mid,10**7))
  if E<0: Lo=mid
  elif E>0:Hi=mid
  else:raise AssertionError((mid,E))
 x=F(Lo,10**7);y=F(Hi,10**7)
 E1,_,_=calc(x,x);E2,_,_=calc(y,y);_,g,_=calc(x,y)
 assert E1<0<E2 and g>0
 print('unique switch kappa in',Lo,'/1e7,',Hi,'/1e7; E endpoints', E1,E2,'derivative',g)
if __name__=='__main__':scan(100,1000)
