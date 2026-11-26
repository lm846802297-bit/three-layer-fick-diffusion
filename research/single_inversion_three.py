"""Interval certificate: exactly one half-spectrum inversion, exactly three reversals."""
import mpmath as m
import numpy as np
from scan_critical import mode_factory
m.mp.dps=76;m.iv.dps=67
M,I=m.mp,m.iv
v,u,z,K='.59796','1.15226','.331107','7.95708'
def layers(A):return [(A.mpf('.5')-A.mpf(z),A.mpf(v)),(A.mpf(z),A.mpf(u)),(A.mpf('.5'),A.mpf(1))]
def calc(q,A,full=False):
    p=A.mpf(K)/q;y=A.mpf(1);norm=0;N=[]
    for j,(w,r) in enumerate(layers(A)):
        a,b,x=y,r*p,q*w*r
        C,S=A.cos(x),A.sin(x)
        if full:
            norm+=w*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*q)
            th=A.atan2(b,a)
            N.append(((x-th-A.pi/2)/A.pi,(-th-A.pi/2)/A.pi))
        y,p=a*C+b*S,(b*C-a*S)/r
        if j==1:mid=y
    return (A.mpf(K)*mid*(1-y)/(norm*q*q),N) if full else p+A.mpf(K)*y/q
def right(q,A):return A.cos(q/2)+A.mpf(K)*A.sin(q/2)/q
def ph_left(q,A):
    (w,v),(z,u),_=layers(A)
    a=q*w*v-A.atan2(v*A.mpf(K),q);c,s=A.cos(a),A.sin(a);r=u/v
    return a+A.atan2((r-1)*c*s,c*c+r*s*s)+q*z*u
def encl(f,guess,eps='1e-48'):
    q=M.findroot(lambda q:f(q,M),guess)
    a,b=[M.nstr(q+d*M.mpf(eps),72) for d in(-1,1)]
    assert f(I.mpf(a),I)*f(I.mpf(b),I)<0
    return I.mpf([a,b])
P=I.mpf
(w,vr),(zz,ur),_=layers(I);kap=P(K)
ellL=w*vr+zz*ur;ellR=P('.5')
assert 59*(ellR-ellL)>1
invs=[]
for n in range(1,30):
    q=encl(right,((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    assert q>(2*n-1)*I.pi and q<2*n*I.pi
    e=ph_left(q,I)-(P(n)-P('.5'))*I.pi
    if e>0:invs.append(n)
    else:assert e<0
assert invs==[2]
print('PASS only half-spectrum inversion: n=2; n>=30 from optical bound',flush=True)
T=P('.0008');qL=ellL/(2*T);qR=ellR/(2*T)
gam=(1+ur/vr)*(1+1/ur)/8; c=kap/(2*gam);g=(ellR**2-ellL**2)/4
d1,d2=w*vr,zz*ur
U=[I.exp(-2*d*qL) for d in (d1,d2,ellR)]
E=(1+U[0])*(1+U[1])*(1+U[2])-1
assert E<1 and (E+U[2])/(1-E)+4*I.sqrt(T)/(ellL*I.sqrt(I.pi))+kap*vr/qL<P('.5')
assert I.exp(-2*ellR*qR)<P(1)/3 and g>T/2
ratio=3*kap*ellL**2*I.sqrt(I.pi)/(c*ellR*I.sqrt(T))*I.exp(-g/T)
assert ratio<P('.56')
print('PASS 0<t<=.0008: right/left < .56',flush=True)
RR=ellR+1/kap; RL=w*vr**2+zz*ur**2+1/kap
F=(RR-RL)/(RR+RL);assert F<0
N=160
_,_,modes=mode_factory(float(u),float(v),float(z),False)
L=[];C=[]
def fl(q):
    a,b=M.mpf(q.a),M.mpf(q.b)
    assert M.floor(a)==M.floor(b) and a>M.floor(a)
    return int(M.floor(a))
for g0 in np.sqrt(modes(float(K),N)[0]):
    g0=M.mpf(float(g0))
    q=encl(calc,(g0-.001,g0+.001)); c,ns=calc(q,I,True)
    assert sum(fl(x)-fl(y) for x,y in ns)==len(L)
    assert not L or q*q>L[-1]
    L.append(q*q);C.append(c)
assert len(L)==N
print('PASS',N,'ordered full-interval modes',flush=True)
# D_min>7/pi^2, kappa<8; |lambda*c| <= 16*(1+2.305*sqrt(lambda)).
# Each omitted term, for t>=.000799, bounded by (100*j)*exp(-7*j*j*t)
# (j >= 160), with geometric ratio <1/2. Also serves for c alone.
assert 1/ur**2>7/I.pi**2 and 2*ur<P('2.305') and kap<8
def tail(t):
    assert t>P('.000798')
    return 200*N*I.exp(-7*N*N*t)
def point(t,d=0):
    t=P(t)
    return (F if d==0 else P(0))-sum(c*(-l)**d*I.exp(-l*t) for l,c in zip(L,C))+P([-1,1])*tail(t)
plan=[('.000799','.028',0,1),('.028','.031',1,-1),('.031','.11',0,-1),
      ('.11','.123',1,1),('.123','1.03',0,1),('1.03','1.13',1,-1),('1.13','1.4',0,-1)]
for aa,bb,d,sg in plan:
    a,b=map(M.mpf,(aa,bb)); W=[-c*(-l)**d for l,c in zip(L,C)]
    stack=[(a,b,0)];passed=0
    while stack:
        a,b,depth=stack.pop()
        mid=(P(a)+P(b))/2; rad=(P(b)-P(a))/2
        Q=[P(0) for _ in range(12)]; Q[0]=F if d==0 else P(0)
        rem=P(0)
        for l,w in zip(L,W):
            f=w*I.exp(-l*mid)
            for j in range(12):
                Q[j]+=f
                f=f*(-l)/(j+1)
            rem+=abs(w)*I.exp(-l*P(a))*(l*rad)**12/M.factorial(12)
        B=sum(abs(Q[i])*rad**i for i in range(1,12))+rem+tail(P(a))
        if sg*Q[0]>B:passed+=1;continue
        assert depth<26,(aa,bb,depth,Q[0],B)
        t=(a+b)/2;stack.extend([(a,t,depth+1),(t,b,depth+1)])
    print('PASS',aa,bb,'derivative' if d else 'value',sg,'cells',passed,flush=True)
for aa,bb,sg in [('.02936168','.02936170',-1),('.11633274','.11633276',1),('1.07769250','1.07769253',-1)]:
    assert sg*point(aa)<0 and sg*point(bb)>0
    print('PASS root bracket',aa,bb,flush=True)
late=sum(abs(c)*I.exp(-l*P('1.4')) for l,c in zip(L,C))+tail(P('1.4'))
assert late<-F
print('PASS t>=1.4 negative. Steady difference',F,flush=True)
