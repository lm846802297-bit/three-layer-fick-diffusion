import numpy as np
from scipy.optimize import brentq
from math import sin,cos
z=.4132;w=.5-z;kappa=6.5;ls=[(w,.04),(z,1.1),(.5,1.)]
Dmin=1/1.1**2
f_inf=-9009/106259009
def calc(k,all=False):
 y=1;p=kappa/k
 norm=0.
 for i,(h,u) in enumerate(ls):
  b=u*p; a=y; x=h*u*k;c=cos(x);s=sin(x)
  if all:norm+=(a*a+b*b)*h/2+(a*a-b*b)*sin(2*x)/(4*u*k)+a*b*sin(x)**2/(u*k)
  y,p=a*c+b*s,(b*c-a*s)/u
  if i==1:at=y
 if all:return kappa*at*(1-y)/norm/k**2
 return p+kappa/k*y
def spectrum(n=100):
 grid=np.arange(.01,4*n+20,.25)
 ks=[]
 for a,b in zip(grid,grid[1:]):
  if calc(a)*calc(b)<0:ks.append(brentq(calc,a,b,xtol=9e-13))
  if len(ks)==n:break
 cs=np.array([calc(k,True) for k in ks])
 return np.array(ks)**2,cs
if __name__=='__main__':
 l,c=spectrum(); print(list(zip(l[:10],c[:10])))
 def val(t,d=0):return (f_inf if d==0 else 0)-np.sum(c*(-l)**d*np.exp(-l*t))
 ts=np.geomspace(.002,5,150)
 for a,b in zip(ts,ts[1:]):
  if val(a)*val(b)<0:
   root=brentq(val,a,b);print('root',root,val(root,1))
 for t in [.002,.003,.004,.005,.006,.007,.008,.009,.01,.02,.1,1,1.1,1.2,1.3,1.4,2]:print(t,val(t),val(t,1))
