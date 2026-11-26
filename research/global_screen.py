from scan import calc
import numpy as np
rng=np.random.default_rng(62026)
best=[]
n=20000
count={}
for i in range(n):
 u=np.exp(rng.uniform(-4,4)); v=np.exp(rng.uniform(-4,4)); k=np.exp(rng.uniform(-6,9))
 if rng.random()<.8:
  u,v=max(u,v,1.0+1e-4),min(u,v,1.0-1e-4)
  lo=.5*(1-v*v)/(u*u-v*v); hi=.5*(1-v)/(u-v)
  if rng.random()<.6: z=lo+(hi-lo)*rng.uniform(-.3,1.4)
  else: z=(lo if rng.random()<.5 else hi)+(hi-lo)*rng.choice([-1,1])*10**rng.uniform(-10,-.2)
 else:z=.5*rng.uniform()
 if not 0<z<.5:continue
 rev=bool(rng.integers(2))
 l,c,t,d=calc(u,v,z,k,n=26,rev=rev)
 o=z*u+(.5-z)*v; R=z*u*u+(.5-z)*v*v
 early=np.sign(.5-o);late=np.sign(.5-R)
 s=np.r_[early,np.sign(d[np.abs(d)>1e-7]),late]
 m=int(np.count_nonzero(s[1:]!=s[:-1]))
 count[m]=count.get(m,0)+1
 if m>=4:
   l,c,t,d=calc(u,v,z,k,n=100,rev=rev); s2=np.r_[early,np.sign(d[np.abs(d)>1e-7]),late];mm=int(np.count_nonzero(s2[1:]!=s2[:-1]))
   if mm>=4:
    sel=np.abs(d)>1e-7; dt=t[sel][np.r_[True,s2[2:-1]!=s2[1:-2]]]; dv=d[sel][np.r_[True,s2[2:-1]!=s2[1:-2]]]
    print('CANDIDATE',i,mm,u,v,z,k,rev,'signs',early,late,list(zip(dt,dv)),flush=True);best.append((u,v,z,k,rev,mm))
 if i%2500==0:print(i,count,flush=True)
print(count,best,flush=True)
