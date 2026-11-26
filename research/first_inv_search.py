from scan_critical import mode_factory
from scipy.optimize import brentq
from math import *
import numpy as np
rng=np.random.default_rng(89326)
def invs(u,v,z,K,rv):
 ls=[(.5-z,v),(z,u)]
 if rv:ls.reverse()
 d=sum(w*a for w,a in ls);gap=.5-d
 if gap<=0:return []
 (w,r),(x,s)=ls
 def phase(q):
  a=q*w*r-atan(K*r/q);c,t=cos(a),sin(a);m=s/r
  return a+atan2((m-1)*c*t,c*c+m*t*t)+q*x*s
 jj=[]
 for n in range(1,ceil((1+1/gap)/2)):
  b=(n-.5)*pi
  q=brentq(lambda q:q/2-atan(K/q)-b,(2*n-1)*pi,2*n*pi)
  if phase(q)>b:jj.append(n)
 return jj
if __name__=='__main__':
 from scan import calc
 ct={}; best=None
 for i in range(70000):
  v=rng.uniform(.13,.93); u=1+exp(rng.uniform(-3.2,.30))
  K=exp(rng.uniform(-.8,5))
  tie=.5*(1-v*v)/(u*u-v*v); opt=.5*(1-v)/(u-v)
  z=tie+10**rng.uniform(-7,-.1)*(opt-tie)
  rv=True if rng.random()<.92 else False # scan: True fast wall
  q=invs(u,v,z,K,not rv)
  if q!=[1]:continue
  L,c,t,d=calc(u,v,z,K,n=32,rev=rv)
  y=np.r_[1,np.sign(d[(abs(d)>3e-7)&(t>.008)]),-1]
  n=int(np.count_nonzero(y[1:]!=y[:-1]))
  ct[n]=ct.get(n,0)+1
  if n>=3:
   rr,_,modes=mode_factory(u,v,z,not rv)
   L,c=modes(K,55)
   R=z*u*u+(.5-z)*v*v+1/K;F=(.5+1/K-R)/(.5+1/K+R)
   d=F-np.exp(-t[:,None]*L)@c
   y=np.r_[1,np.sign(d[(abs(d)>1e-10)&(t>.006)]),-1]
   n=int(np.count_nonzero(y[1:]!=y[:-1]))
   if n>=3:
    print('FOUND',u,v,z,K,rv,q,n,flush=True)
    print([(ti,di) for ti,di in zip(t,d) if ti>.006 and abs(di)>1e-10][::10],flush=True)
    break
  if i%4000==0:print(i,ct,flush=True)
 print(i,ct,flush=True)
