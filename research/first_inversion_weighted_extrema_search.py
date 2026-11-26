"""Floating-point search for possible finite-time tangencies of W."""
from scan_critical import mode_factory
from one_inv_search import invs
from scipy.optimize import brentq
import numpy as np
from math import exp
rng=np.random.default_rng(511362)
t=np.geomspace(.0002,1.5,190)
def extrema(u,v,z,K,rv,N):
    L,C=mode_factory(u,v,z,rv)[2](K,N)
    G=L[1:]-L[1]; W=(C*L*(L[0]-L))[1:]
    Wd=-G*W
    ts=t[t*G[-1]>48]
    E=np.exp(-ts[:,None]*G); d=E@Wd
    size=E@abs(Wd)
    good=(size>1e-7)&(abs(d)>1e-7*size)
    ts,d=ts[good],d[good]
    results=[]
    for i in np.flatnonzero(d[:-1]*d[1:]<0):
       T=brentq(lambda x:np.exp(-x*G)@Wd,ts[i],ts[i+1],xtol=1e-14)
       E=np.exp(-T*G); mass=E@abs(W); y=E@W
       if mass>1e-6:
           results.append((float(T),float(y/mass),bool(d[i]<0),float(y)))
    return results
if __name__=='__main__':
    counts={};best={}; samples=0
    for i in range(18500):
        u=1+exp(rng.uniform(-4,4));v=exp(rng.uniform(-8,-.007));K=exp(rng.uniform(-4,8))
        tie=.5*(1-v*v)/(u*u-v*v);z=tie*(1-10**rng.uniform(-9,-.015))
        rv=bool(rng.integers(2))
        if .5-(.5-z)*v-z*u<.005 or invs(u,v,z,K,rv)!=[1]:continue
        samples+=1
        p=(u,v,z,K,rv)
        ex=extrema(*p,85)
        pat=tuple(('min' if a else 'max','+' if r>0 else '-') for _,r,a,_ in ex)
        counts[pat]=counts.get(pat,0)+1
        for T,r,a,y in ex:
            if abs(r)<abs(best.get(a,(0,1,0))[1]):
                best[a]=(T,r,y,p)
        if any(a and r>0 for _,r,a,_ in ex):
            print('POSITIVE MINIMUM',p,extrema(*p,145),flush=True)
            break
    print('samples',samples,'patterns',counts,'nearest',best,flush=True)
