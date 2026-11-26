"""Float-only targeted screen for >=4 reversals: first-only inversion and d>0."""
from one_inv_search import invs
from scipy.optimize import brentq
from scan_critical import mode_factory
from scan import calc
import numpy as np
rng=np.random.default_rng(94372)
counts={}; found=[]
for i in range(35000):
    u=1+np.exp(rng.uniform(-4,3.5));v=np.exp(rng.uniform(-6,-.03));K=np.exp(rng.uniform(-3.5,7))
    tie=.5*(1-v*v)/(u*u-v*v)
    z=tie*(1-10**rng.uniform(-9,-.01))
    # rev=True means slow layer at wall in the Prüfer helper.
    rev=bool(rng.integers(2))
    if invs(u,v,z,K,rev)!=[1]:continue
    _,_,t,d=calc(u,v,z,K,30,not rev)
    start=.008
    part=(t>=start)&(abs(d)>1e-7)
    signs=np.r_[1,np.sign(d[part]),1]
    n=int(np.count_nonzero(signs[:-1]!=signs[1:]))
    counts[n]=counts.get(n,0)+1
    if n<4:continue
    # Use exact piecewise eigenfunction formula with floating-point roots.
    _,_,modes=mode_factory(u,v,z,rev)
    L,C=modes(K,80)
    RR=.5+1/K;RL=z*u*u+(.5-z)*v*v+1/K
    D=(RR-RL)/(RR+RL)-np.exp(-t[:,None]*L)@C
    active=(t>=start)&(abs(D)>1e-9)&(t*L[-1]>28)
    sig=np.r_[1,np.sign(D[active]),1]
    if np.count_nonzero(sig[:-1]!=sig[1:])>=4:
       found.append((u,v,z,K,rev));print('candidate',found[-1],flush=True)
print('screen approximate counts',counts,'exact-transfer candidates',found)
