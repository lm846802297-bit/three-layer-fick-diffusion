from first_inversion_weighted_search import run
from one_inv_search import invs
import numpy as np
from math import exp
rng=np.random.default_rng(16592)
ok=0
for i in range(20000):
    rr=10**rng.uniform(-4,2)
    v=1/(1+rr) if i%2==0 else exp(rng.uniform(-12,-.001))
    u=1+rr*10**rng.uniform(-2,2)
    K=10**rng.uniform(-3,5)
    tie=.5*(1-v*v)/(u*u-v*v)
    z=tie*(1-10**rng.uniform(-13,-.015))
    rv=bool(rng.integers(2))
    if .5-(.5-z)*v-z*u<.005:continue # cap the finite half-spectrum screen
    if invs(u,v,z,K,rv)!=[1]:continue
    ok+=1
    rep=run(u,v,z,K,rv,75)
    if sum(w<0 for _,w in rep)>=2:
        rep=run(u,v,z,K,rv,155)
        if sum(w<0 for _,w in rep)>=2:
            print('FOUND',i,ok,u,v,z,K,rv,rep,flush=True)
            break
print('accepted',ok)
