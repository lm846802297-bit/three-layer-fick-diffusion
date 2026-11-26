"""Uniform interval certificate for an isolated finite-time fold slice."""
import mpmath as m
import numpy as np
from scan_critical import mode_factory
M,I=m.mp,m.iv
M.dps=73;I.dps=65
P=I.mpf
ka,kb='5.0902282130','5.0902282140'
def show(x):return M.nstr(x,70)
def calc(q,K,A,full=False):
    M0=A.mpf
    p=K/q;y=M0(1);norm=M0(0);N=[]
    for j,(ww,rr) in enumerate((('.38','.02'),('.12','2'),('.5','1'))):
        w,r=M0(ww),M0(rr)
        a,b,x=y,r*p,q*w*r; C,S=A.cos(x),A.sin(x)
        if full:
            norm+=w*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*q)
            th=A.atan2(b,a);N.append(((x-th-A.pi/2)/A.pi,(-th-A.pi/2)/A.pi))
        y,p=a*C+b*S,(b*C-a*S)/r
        if j==1:mid=y
    return (K*mid*(1-y)/(norm*q*q),N) if full else p+K*y/q
def eigbr(guess,K,eps='1e-47'):
    q=M.findroot(lambda q:calc(q,K,M),(guess-.001,guess+.001))
    a,b=[show(q+d*M.mpf(eps)) for d in (-1,1)]
    assert calc(P(a),P(show(K)),I)*calc(P(b),P(show(K)),I)<0
    return a,b
def cocoe(aa,bb,N):
    a,b=map(M.mpf,(aa,bb));K=P([aa,bb]); mid=(a+b)/2
    _,_,modes=mode_factory(2,.02,.12,False)
    L=[];C=[]
    def fl(q):
        u,v=M.mpf(q.a),M.mpf(q.b);assert M.floor(u)==M.floor(v) and u>M.floor(u)
        return int(M.floor(u))
    for g0 in np.sqrt(modes(float(mid),N)[0]):
        g=M.mpf(float(g0));lo,e=eigbr(g,a);f,hi=eigbr(g,b)
        lo,hi=P(lo),P(hi);q=P([M.mpf(lo.a),M.mpf(hi.b)])
        c,ns=calc(q,K,I,True)
        assert sum(fl(x)-fl(y) for x,y in ns)==len(L)
        assert not L or q*q>L[-1]
        L.append(q*q);C.append(c)
    assert len(L)==N
    return L,C,K
ell1=P('.38')*P('.02');ell2=P('.12')*2;ell=ell1+ell2
RR=P('.5');rL=P('.38')*P('.02')**2+P('.12')*4
K=P([ka,kb]);gap=RR-ell
assert 5*gap>1
def right(q,k,A):return A.cos(q/2)+k*A.sin(q/2)/q
sl=[]
for n in [1,2]:
    q=M.findroot(lambda x:right(x,(M.mpf(ka)+M.mpf(kb))/2,M),((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    a,b=[P(show(q+d*M.mpf('1e-8'))) for d in (-1,1)]
    assert a>(2*n-1)*I.pi and b<2*n*I.pi
    assert right(a,K,I)*right(b,K,I)<0
    q=P([M.mpf(a.a),M.mpf(b.b)])
    A=ell1*q-I.atan2(K*P('.02'),q);c,t=I.cos(A),I.sin(A)
    g=A+I.atan2(99*c*t,c*c+100*t*t)+ell2*q-(n-.5)*I.pi
    sl.append(g)
assert sl[0]>0 and sl[1]<0
print('PASS only the first half-spectrum mode is inverted, uniformly in kappa bracket',flush=True)
T=P('.0018'); q=ell/(2*T);qr=RR/(2*T)
gm=(1+100)*(1+P('.5'))/8; C=K/(2*gm)
E=[I.exp(-2*d*q) for d in (ell1,ell2,RR)]
et=(1+E[0])*(1+E[1])*(1+E[2])-1
g=(RR**2-ell**2)/4
assert (et+E[2])/(1-et)+2/(q*I.sqrt(I.pi*T))+K*P('.02')/q<P('.95')
assert I.exp(-2*RR*qr)<P(1)/3 and g>T/2
# A_L lower bound is C/(20*q^2*sqrt(pi*t))*exp(-ell^2/(4t)).
# A_R upper bound is 3*kappa/q_R*exp(-.5^2/(4t)).
assert 30*K*ell**2*I.sqrt(I.pi)/(C*RR*I.sqrt(T))*I.exp(-g/T)<1
print('PASS uniform 0<t<=.0018 positive',flush=True)
def tod(K):
    return (RR-rL)/(RR+rL+2/K)
def tail(t,N,d):
    assert t>P('.001798')
    # lambda_j>=2*j^2, |lambda*c| <= 12*(1+4 sqrt(lambda)).
    # Each n>N is bounded at lower eigenvalue; ratio <1/2.
    return (160*N if d<2 else 400*N**3)*I.exp(-2*N*N*t)
def point(t,L,C,K,d=0):
    t=P(t)
    return (tod(K) if d==0 else P(0))-sum(c*(-l)**d*I.exp(-l*t) for l,c in zip(L,C))+P([-1,1])*tail(t,len(L),d)
L,C,K=cocoe(ka,kb,160)
print('PASS 160 indexed full-spectrum modes across parameter bracket',flush=True)
def slab(aa,bb,d):
    W=[-c*(-l)**d for l,c in zip(L,C)]
    stk=[(M.mpf(aa),M.mpf(bb),0)];cnt=0
    while stk:
        a,b,depth=stk.pop();mid=(P(a)+P(b))/2; rad=(P(b)-P(a))/2
        Q=[P(0)]*12;Q[0]=tod(K) if d==0 else P(0);rem=P(0)
        for l,w in zip(L,W):
            f=w*I.exp(-l*mid)
            for j in range(12):
                Q[j]+=f; f=f*(-l)/(j+1)
            rem+=abs(w)*I.exp(-l*P(a))*(l*rad)**12/M.factorial(12)
        B=sum(abs(Q[j])*rad**j for j in range(1,12))+rem+tail(P(a),len(L),d)
        if Q[0]>B:cnt+=1;continue
        assert depth<26,(aa,bb,a,b,Q[0],B)
        x=(a+b)/2;stk.extend([(a,x,depth+1),(x,b,depth+1)])
    print('PASS time slab',aa,bb,'order',d,'positive on',cnt,'pieces',flush=True)
for a,b,d in [('.001799','.1519',0),('.1519','.1526',2),('.1526','1.2',0)]:slab(a,b,d)
assert point('.1519',L,C,K,1)<0 and point('.1526',L,C,K,1)>0
assert sum(abs(c)*I.exp(-l*P('1.2')) for l,c in zip(L,C))+tail(P('1.2'),len(L),0)<tod(K)
print('PASS exterior times positive and unique minimum inside time window',flush=True)
L,C,K=cocoe(ka,ka,18); x=point('.1522445625',L,C,K)
assert x<0
print('PASS kappa interval lower end has two strict crossings',x,flush=True)
L,C,K=cocoe(kb,kb,18); x=point('.1522445625',L,C,K)
assert point('.15224455',L,C,K,1)<0 and point('.15224458',L,C,K,1)>0
# Use a derivative bound on the narrow interval from convexity.
bound=max(abs(point('.15224455',L,C,K,1)),abs(point('.15224458',L,C,K,1)))*P('3e-8')
assert x>bound
print('PASS kappa interval upper end has no crossings',x,bound,flush=True)
