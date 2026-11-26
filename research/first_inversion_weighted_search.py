from scan_critical import mode_factory
from one_inv_search import invs
import numpy as np
from math import exp
rng=np.random.default_rng(8794181)
t=np.geomspace(.000035,2,460)
def run(u,v,z,K,rv,N):
    _,_,md=mode_factory(u,v,z,rv)
    L,C=md(K,N)
    coef=C*L*(L[0]-L)
    v=np.exp(-t[:,None]*(L-L[1]))@coef
    good=(t*(L[-1]-L[1])>43)&(abs(v)>max(1e-7,1e-4*abs(coef[1])))
    xs=np.flatnonzero(good)
    xs=xs[np.r_[True,np.diff(np.sign(v[xs]))!=0]] if len(xs) else xs
    return [(float(t[i]),float(v[i])) for i in xs]
if __name__=='__main__':
    ok=0
    for i in range(22000):
        u=1+exp(rng.uniform(-3,3.5));v=exp(rng.uniform(-7,-.018));K=exp(rng.uniform(-3,7))
        tie=.5*(1-v*v)/(u*u-v*v)
        z=tie*(1-10**rng.uniform(-6,-.06))
        rv=bool(rng.integers(2))
        if invs(u,v,z,K,rv)!=[1]:continue
        ok+=1
        rep=run(u,v,z,K,rv,70)
        if sum(w<0 for _,w in rep)>=2:
            rep=run(u,v,z,K,rv,155)
            if sum(w<0 for _,w in rep)>=2:
                print('FOUND',i,ok,u,v,z,K,rv,rep,flush=True)
                break
    print('accepted',ok)
