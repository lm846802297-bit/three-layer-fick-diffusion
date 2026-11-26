from scan import calc
import numpy as np
from scipy.optimize import differential_evolution
def objective(x,n=35,out=False):
 u,v,k,f=x
 lb=(1-v*v)/(u*u-v*v); ub=(1-v)/(u-v)
 z=.5*(lb+f*(ub-lb)); L,c,t,d=calc(u,v,z,k,n,rev=True)
 ix=np.argmax(np.minimum(-np.minimum.accumulate(d),d))
 a=np.argmin(d[:ix+1])
 score=min(-d[a],d[ix])
 return (score,z,t[a],t[ix],d[a],d[ix]) if out else -score
if __name__=='__main__':
 r=differential_evolution(objective,[(1.003,1.6),(.002,.2),(.2,35),(0,.6)],seed=11,popsize=12,maxiter=55,tol=.001)
 print(r.x,r.fun); print(objective(r.x,200,True))
 for x in [[1.1,.04,6.5,.05],[1.1,.02,5,.12],[1.12,.02,6,.05], [1.08,.025,5,.1]]:
  print(x,objective(x,150,True))
