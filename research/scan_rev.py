import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
# mesh piecewise, cell centers; Robin with series half-cell resistance
def calc(u,v,z,k,n=50,rev=False):
 hs=np.repeat(np.array([z,.5-z,.5])/n,n); ds=np.repeat([u**-2,v**-2,1.],n)
 if rev: hs[:2*n]=np.r_[hs[n:2*n],hs[:n]];ds[:2*n]=np.r_[ds[n:2*n],ds[:n]]
 g=1/(hs[:-1]/(2*ds[:-1])+hs[1:]/(2*ds[1:]))
 gl,gr=1/(np.array([hs[0]/(2*ds[0]),hs[-1]/(2*ds[-1])])+1/k)
 h=np.r_[gl,g]+np.r_[g,gr]
 vals,V=eigh_tridiagonal(h/hs,-g/np.sqrt(hs[:-1]*hs[1:]))
 src=(V[2*n-1]/np.sqrt(hs[2*n-1])+V[2*n]/np.sqrt(hs[2*n]))/2
 c=(V[0]*gl/np.sqrt(hs[0])-V[-1]*gr/np.sqrt(hs[-1]))*src/vals
 # analytic steady removed from numerical cancellation; expression for short time near zero only approximate
 Rl=z*u*u+(.5-z)*v*v+1/k; Rr=.5+1/k
 F=(Rr-Rl)/(Rr+Rl)
 ts=np.geomspace(0.001,max(20,15/vals[0]),255)
 delta=F-np.exp(-ts[:,None]*vals)@c
 return vals,c,ts,delta
if __name__=='__main__':
 rng=np.random.default_rng(561)
 best=[]; stats={}
 for i in range(13000):
  u=1+10**rng.uniform(-1,2.2); v=10**rng.uniform(-1.7,-.0001);k=10**rng.uniform(-2.5,4)
  l=(1-v*v)/(u*u-v*v); h=(1-v)/(u-v); f=rng.choice([rng.uniform(),10**rng.uniform(-8,0),1-10**rng.uniform(-8,0)])
  z=.5*(l+f*(h-l)); lam,c,t,d=calc(u,v,z,k,n=30,rev=True)
  signs=np.sign(d[abs(d)>1e-6]); ups=np.sum(signs[1:]>signs[:-1]); downs=np.sum(signs[1:]<signs[:-1]);
  st=(int(ups),int(downs)); stats[st]=stats.get(st,0)+1
  if ups and downs>=1:
   ll,cc,tt,dd=calc(u,v,z,k,n=120,rev=True); mask=abs(dd)>1e-7; s=np.sign(dd[mask]); 
   if np.sum(s[1:]!=s[:-1])>=3:
    print('HIT', u,v,z,k,[(float(x),float(y)) for x,y in zip(tt[mask][np.r_[True,s[1:]!=s[:-1]]], dd[mask][np.r_[True,s[1:]!=s[:-1]]])],flush=True); break
 print(i,stats)
