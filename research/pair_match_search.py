from math import *
from scipy.optimize import brentq
from numpy.random import default_rng
rng=default_rng(514)
for i in range(50000):
 u=1+exp(rng.uniform(-4,1));r=exp(rng.uniform(-3,-.01));K=exp(rng.uniform(-3,7))
 x=.5*(1-r*r)/(u*u-r*r);y=.5*(1-r)/(u-r)
 z=x+(y-x)*rng.uniform(.005,.85)
 gap=.5-(.5-z)*r-z*u
 if gap<.025:continue
 rv=rng.random()<.25
 w,a,ww,b=(z,u,.5-z,r) if rv else (.5-z,r,z,u)
 def ph(k):
  q=k*w*a-atan(K*a/k);c,s=cos(q),sin(q);m=b/a
  return q+atan2((m-1)*c*s,c*c+m*s*s)+k*ww*b
 L=[];R=[]
 for n in range(1,ceil((1+1/gap)/2)):
  g=(n-.5)*pi
  rr=brentq(lambda q:q/2-atan(K/q)-g,(2*n-1)*pi,2*n*pi)
  ll=brentq(lambda q:ph(q)-g,.001,max(40,rr*4))
  L.append(ll*ll);R.append(rr*rr)
 j=[i for i,(a,b) in enumerate(zip(L,R)) if a<b]
 if len(j)==1 and j[0]>0:
  i=j[0]; choices=[(L[k]-R[k]-(R[i]-L[i]),k) for k in range(i)]
  margin,k=max(choices)
  if margin>0.1:
    print('FOUND',r,u,z,K,rv,'inverted',i+1,'partner',k+1,'margin',margin,L[:i+1],R[:i+1],flush=True)
    break
else:print('NOT FOUND')
