import numpy as np
from scipy.optimize import brentq
from scan_critical import mode_factory
from one_inv_search import invs
rr=0; good=0
for v in np.geomspace(.02,.94,50):
 for u in 1+np.geomspace(.027,5,45):
    z=.5*(1-v*v)/(u*u-v*v)
    for rev in (False,True):
      root,co,modes=mode_factory(u,v,z,rev)
      def f(k):return co(root(0,k),k)
      ks=np.geomspace(.001,1e4,16)
      for a,b in zip(ks,ks[1:]):
        if f(a)*f(b)>=0:continue
        K=brentq(f,a,b); rr+=1
        ii=invs(u,v,z,K,rev)
        if ii==[1]:
          good+=1
          print('found',u,v,z,K,rev,ii,flush=True)
          L,C=modes(K,60)
          t=np.geomspace(.003,.65,220);d=-np.exp(-t[:,None]*L[1:])@C[1:]
          print('c2',C[1],'min',np.min(d),'sign witnesses',[(s,np.interp(s,t,d)) for s in [.006,.02,.05,.1,.2,.4]],flush=True)
          if np.min(d)<-1e-8:raise SystemExit
          if good>=15:raise SystemExit
print('totals',rr,good)
