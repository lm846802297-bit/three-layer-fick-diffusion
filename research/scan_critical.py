"""Exploratory float screen only, NOT a proof or interval calculation."""
import numpy as np
from scipy.optimize import brentq
from math import sin,cos,atan2,atan,pi
def mode_factory(u,v,z,rev):
 ls=[(.5-z,v),(z,u),(.5,1.)]
 if rev: ls[:2]=ls[:2][::-1]
 def phase(k,K):
    ang=-atan(ls[0][1]*K/k)
    for i,(h,r) in enumerate(ls):
        ang+=h*r*k
        if i<2:
            x,y=cos(ang),sin(ang)
            yy=y*ls[i+1][1]/r
            ang+=atan2(x*(yy-y),x*x+y*yy)
    return ang-atan(K/k)
 def root(n,K):
    hi=4*(n+1)
    while phase(hi,K)<n*pi:hi*=2
    return brentq(lambda k:phase(k,K)-n*pi,1e-9,hi,xtol=1e-14)
 def coeff0(k,K):
    y,p=1,K/k
    for h,r in ls:
        a,b,x=y,r*p,h*r*k
        y,p=a*cos(x)+b*sin(x),(b*cos(x)-a*sin(x))/r
    return 1-y
 def modes(K,N):
    l=[];c=[]
    for n in range(N):
        k=root(n,K); y,p=1,K/k;norm=0
        for i,(h,r) in enumerate(ls):
            a,b,x=y,r*p,h*r*k
            co,si=cos(x),sin(x)
            norm+=(a*a+b*b)*h/2+(a*a-b*b)*sin(2*x)/(4*r*k)+a*b*si*si/(r*k)
            y,p=a*co+b*si,(b*co-a*si)/r
            if i==1:mid=y
        l.append(k*k);c.append(K*mid*(1-y)/norm/k**2)
    return np.array(l),np.array(c)
 return root,coeff0,modes

def main():
    counts={}; nc=0
    for rev in [False,True]:
        for u in 1+np.geomspace(.025,8,18):
            for v in np.unique(np.r_[np.geomspace(.012,.94,16),np.linspace(.15,.97,8)]):
                z=.5*(1-v*v)/(u*u-v*v)
                root,co,modes=mode_factory(u,v,z,rev)
                def f(K):return co(root(0,K),K)
                ks=np.geomspace(.001,1e4,21)
                ss=[f(k) for k in ks]
                for j in range(len(ks)-1):
                    if not ss[j]*ss[j+1]<0:continue
                    nc+=1
                    K=brentq(f,ks[j],ks[j+1],xtol=1e-13)
                    L,C=modes(K,N=55)
                    if C[1]<=1e-7:continue
                    t=np.geomspace(.003,1.5,260)
                    scaled=-np.exp(-t[:,None]*(L[None,1:]-L[1]))@C[1:]
                    s=np.r_[1,np.sign(scaled[abs(scaled)*np.exp(-L[1]*t)>1e-9]),-1]
                    crossings=int(np.sum(s[1:]!=s[:-1]))
                    counts[(rev,crossings)]=counts.get((rev,crossings),0)+1
                    if crossings>1: print('needs validation',u,v,z,K,rev,crossings,flush=True)
    print('Candidate roots',nc,'with c2>1e-7 classified by (reverse order, early sign changes):',counts)
if __name__=='__main__':main()
