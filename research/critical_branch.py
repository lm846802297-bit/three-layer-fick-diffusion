"""Reduced equations and exploratory continuation."""
import numpy as np
from scipy.optimize import root
from math import sin,cos,sqrt,pi,tan
from scan_critical import mode_factory

def eq(v,r,rev=False):
 h,k1,k2=v
 # analytic parameterization, with no use outside the physical domain
 if not 0<h<.5 or not 0<k1<pi or not pi<k2<2*pi:
    return [1e3+abs(h),1e3+abs(k1),1e3+abs(k2)]
 z=.5-h; R=.5-h*r*r; u=sqrt(R/z); K=-k2/tan(k2/2)
 def l(k):
    a,b=1,K/k
    seq=[(h,r),(z,u)]
    for w,s in seq[::-1] if rev else seq:
       x=w*s*k; c=cos(x); t=sin(x)
       a,b=c*a+s*t*b,-t*a/s+c*b
    return a,b
 a,b=l(k1); c=cos(k1/2)+K/k1*sin(k1/2)
 d=-sin(k1/2)+K/k1*cos(k1/2)
 return [a*d+b*c,a-c,l(k2)[0]]

def cont():
 from scipy.optimize import root
 mid=[.20448868701134,sqrt(7.18405163448),sqrt(26.03766133788)]
 vals=[]; misses=[]
 for seq in [np.geomspace(.5,.001,110),np.linspace(.5,.995,105)]:
    x=mid
    for r in seq:
       rt=root(eq,x,args=(r,),tol=1e-10)
       if np.max(np.abs(rt.fun))>1e-9:
          misses.append(r);continue
       x=rt.x
       h,k1,k2=x;z=.5-h; u=sqrt((.5-h*r*r)/z); K=-k2/tan(k2/2)
       fr,_,modes=mode_factory(u,r,z,False)
       L,C=modes(K,65)
       if max(abs(C[:2]))>1e-7 or abs(L[1]-k2*k2)>1e-6:raise RuntimeError('misindexed')
       t=np.geomspace(.005,1.5,380)
       y=-np.exp(-t[:,None]*(L[None,2:]-L[2]))@C[2:]
       vals.append([r,h,k1,k2,K,C[2],np.min(y),t[np.argmin(y)]])
 a=np.array(vals)
 np.savetxt('critical_branch_output.txt',a,header='r h k1 k2 kappa c3 min_scaled_delta t_at_min (floating-point search only)')
 print('successful',len(vals),'misses',misses)
 print('minmax c3',min(a[:,5]),max(a[:,5]),'min scaled delta',min(a[:,6]))
 for r in [.001,.01,.05,.1,.2,.5,.8,.95,.995]:
   print(a[np.argmin(abs(a[:,0]-r))])
if __name__=='__main__':cont()
