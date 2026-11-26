from scipy.optimize import brentq,minimize_scalar
from numpy.random import default_rng
from math import *
rng=default_rng(6072)
for j in range(70000):
 u=1+exp(rng.uniform(-4,1.3));v=exp(rng.uniform(-3,-.015));k=exp(rng.uniform(-3.5,6))
 tie=.5*(1-v*v)/(u*u-v*v); opt=.5*(1-v)/(u-v)
 z=tie+(opt-tie)*rng.uniform(.001,.75);w=.5-z;gap=.5-(w*v+z*u)
 if gap<.035: continue
 seq=[(w,v),(z,u)]
 rev=rng.random()<.25
 if rev:seq.reverse()
 (w1,r1),(w2,r2)=seq
 def ph(q):
  a=q*w1*r1-atan(k*r1/q);x,y=cos(a),sin(a);r=r2/r1
  return a+atan2((r-1)*x*y,x*x+r*y*y)+q*w2*r2
 L=[];R=[]
 N=max(5,ceil((1+1/gap)/2)-1)
 for n in range(1,N+1):
  q=brentq(lambda x:x/2-atan(k/x)-(n-.5)*pi,(2*n-1)*pi,2*n*pi)
  l=brentq(lambda x:ph(x)-(n-.5)*pi,1e-8,max(10,6*q))
  L.append(l*l);R.append(q*q)
 inv=[i for i in range(N) if L[i]<R[i]]
 if len(inv)!=1:continue
 i=inv[0];a,b=L[i],R[i]
 lo=sum(L[h]-R[h] for h in range(i))
 if i==0 or lo>=.97*(b-a):continue
 def rat(t):
  if t==0:return sum(L[h]-R[h] for h in range(N) if h!=i)/(b-a)
  return sum(exp((a-R[h])*t)*(-expm1(-(L[h]-R[h])*t)) for h in range(N) if h!=i)/(-expm1(-(b-a)*t))
 ans=minimize_scalar(rat,bounds=(0,min(1,20/(a-R[0]))),method='bounded')
 if ans.fun>1.08:
  print(v,u,z,k,rev,'inverse',i+1,'N',N,'lo/width',lo/(b-a),'lower min',ans.fun,ans.x, 'spec',list(zip(L,R)),flush=True)
  break
else: print('none')
