import numpy as np
from scipy.optimize import brentq,root
from math import sin,cos
def make(kappa,z=.4132,u=1.1,v=.04,n=32,rev=False):
 ls=[(.5-z,v),(z,u),(.5,1.)]
 if rev: ls[:2]=ls[:2][::-1]
 rr=.5+1/kappa; rl=z*u*u+(.5-z)*v*v+1/kappa; F=(rr-rl)/(rr+rl)
 def f(k,full=False):
  y,p=1,kappa/k; norm=0
  for i,(h,u) in enumerate(ls):
   b=u*p;a=y; x=h*u*k;c=cos(x);s=sin(x)
   if full:norm+=(a*a+b*b)*h/2+(a*a-b*b)*sin(2*x)/(4*u*k)+a*b*s*s/(u*k)
   y,p=a*c+b*s,(b*c-a*s)/u
   if i==1:at=y
  return kappa*at*(1-y)/norm/k**2 if full else p+kappa/k*y
 grid=np.arange(.0001,4*n+20,.15); ks=[]
 for a,b in zip(grid,grid[1:]):
  if f(a)*f(b)<0:ks.append(brentq(f,a,b,xtol=1e-14))
  if len(ks)==n:break
 l=np.array(ks)**2;c=np.array([f(k,True) for k in ks])
 def val(t,d=0):return (F if d==0 else 0)-np.sum(c*(-l)**d*np.exp(-l*t))
 return val,l,c
if __name__=='__main__':
 for k in [.01,.1,.3,1,2,3,4,5,6,6.5,7,8,10,15,20,50,100]:
  val,l,c=make(k); T=np.geomspace(.004,10,120); print(k,[(round(brentq(val,a,b),10)) for a,b in zip(T,T[1:]) if val(a)*val(b)<0], c[0])
 for x in [[.1,4],[.05,20],[.02,2],[.2,2],[.025,20]]:
  def obj(a):
   t,k=np.exp(a); f,*_=make(k,n=20); return [1000*f(t),100*f(t,1)]
  r=root(obj,np.log(x)); print(r.success,np.exp(r.x),r.fun)
 z=.5*(1-.04**2)/(1.1**2-.04**2)
 for k in [.01,.1,1,6.5,10,100,1000]:
  print('tie',k,[(make(k,z=z,rev=rev)[2][0]) for rev in [False,True]])
