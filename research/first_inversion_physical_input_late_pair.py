"""Interval proof of the critical point d=c2=0 and strict capture positivity."""
import mpmath as mp
from late_triple import AD
from scan_critical import mode_factory
M,I=mp.mp,mp.iv
M.dps=80; I.dps=74
P=I.mpf
us,vs='1.8873913557547013','0.12633315782720467'
def state(q,z,k,A,kind,ad=False):
    f=AD if ad else A.mpf
    u,v=f(us),f(vs); half=f('.5')
    layers=[(half-z,v),(z,u),(half,f(1))]
    if kind=='right':layers=layers[2:]
    if kind=='left':layers=layers[:2]
    y,p=f(1),k/q
    J=[];norm=f(0)
    for j,(h,r) in enumerate(layers):
        a,b,x=y,r*p,q*h*r
        S,C=(x.sin(),x.cos()) if ad else (A.sin(x),A.cos(x))
        if kind=='modal':
            g=A.atan2(b,a)
            J.append(((x-g-A.pi/2)/A.pi,(-g-A.pi/2)/A.pi))
            norm+=h*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*q)
        y,p=a*C+b*S,(b*C-a*S)/r
        if j==1:mid=y
    if kind=='modal':return y,p,mid,J,k*mid*(1-y)/(norm*q*q)
    return y,p
def eq(args,A):
    z,k,a,b=args
    ad=isinstance(z,AD)
    l=state(a,z,k,A,'left',ad)[0]; r=state(a,z,k,A,'right',ad)[0]
    y,p=state(b,z,k,A,'full',ad)
    f=AD if ad else A.mpf
    u,v=f(us),f(vs)
    return l,r,p+k*y/b,(f('.5')-z)*v*v+z*u*u-f('.5')
def s(x):return M.nstr(x,77)
X0=M.findroot(lambda *xs:eq(xs,M),('.13874232306','9.75011478186','5.289104','7.875016'))
eps=M.mpf('1e-54')
box=[P([s(x-eps),s(x+eps)]) for x in X0]; mid=[P(s(x)) for x in X0]
Y0=M.inverse(M.matrix(M.jacobian(lambda *xs:eq(xs,M),X0)))
Y=[[P(s(Y0[i,j])) for j in range(4)] for i in range(4)]
J=[r.d for r in eq([AD(x,[P(int(i==j)) for j in range(4)]) for i,x in enumerate(box)],I)]
F=eq(mid,I)
for i in range(4):
    R=[int(i==j)-sum(Y[i][l]*J[l][j] for l in range(4)) for j in range(4)]
    assert sum(abs(x) for x in R)<P('.001')
    assert abs(-sum(Y[i][j]*F[j] for j in range(4))+sum(R[j]*(box[j]-mid[j]) for j in range(4)))<P(s(eps))/2
print('PASS existence and full rank',*(M.nstr(x,24) for x in X0))
z,k,a,b=box; u,v=P(us),P(vs);h=P('.5')-z
assert I.pi<a<2*I.pi
yb,_=state(a,z,k,I,'left')
t=a*h*v-I.atan2(k*v,a); c,d=I.cos(t),I.sin(t)
def phase(q):
    t=q*h*v-I.atan2(k*v,q);c,d=I.cos(t),I.sin(t)
    return t+I.atan2((u/v-1)*c*d,c*c+(u/v)*d*d)+q*z*u
assert I.pi/4<phase(a)<3*I.pi/4
for n in range(2,50):
    gap=P('.5')-h*v-z*u
    if (2*n-1)*gap>=1:break
    q=M.findroot(lambda q:q/2-M.atan2(X0[1],q)-(n-M.mpf('.5'))*M.pi,((2*n-1)*M.pi,2*n*M.pi))
    Q=P([s(q-M.mpf('1e-50')),s(q+M.mpf('1e-50'))])
    assert (2*n-1)*I.pi<Q<2*n*I.pi
    p1=P(M.mpf(Q.a));p2=P(M.mpf(Q.b))
    assert (p1/2-I.atan2(k,p1)-(n-P('.5'))*I.pi)*(p2/2-I.atan2(k,p2)-(n-P('.5'))*I.pi)<0
    assert phase(Q)<(n-P('.5'))*I.pi
else:raise AssertionError()
assert (h*v*v+z*u*u-P('.5')).a<=0<=(h*v*v+z*u*u-P('.5')).b
print('PASS first half-spectrum tie, higher modes strict, static tie')
def floor(Q):
    lo,hi=M.mpf(Q.a),M.mpf(Q.b)
    assert M.floor(lo)==M.floor(hi) and lo>M.floor(lo)
    return int(M.floor(lo))
def index(J):return sum(floor(x)-floor(y) for x,y in J)

_,_,X,J,c3=state(b,z,k,I,'modal')
assert X<0 and index(J)==2 and c3>0
def full_r(q,z,k,A):
    y,p=state(q,z,k,A,'full');return p+k*y/q
def bracket(fn,guess):
    q=M.findroot(fn,(M.mpf(guess)-M.mpf('.01'),M.mpf(guess)+M.mpf('.01')))
    aa,bb=[P(s(q+d*M.mpf('1e-50'))) for d in (-1,1)]
    return aa,bb,P([M.mpf(aa.a),M.mpf(bb.b)])
g,j,Q=bracket(lambda q:full_r(q,*X0[:2],M),'3.085467')
assert full_r(g,z,k,I)*full_r(j,z,k,I)<0
_,_,_,J,c1=state(Q,z,k,I,'modal')
assert index(J)==0 and c1<0
g,j,U=bracket(lambda q:state(q,*X0[:2],M,'right')[0],'10.0')
assert state(g,z,k,I,'right')[0]*state(j,z,k,I,'right')[0]<0
assert 3*I.pi<U<4*I.pi
lam=Q**2; m=U**2; br=a**2
R=state(Q,z,k,I,'left')[0]/state(Q,z,k,I,'right')[0]
hi=br*(m-lam)/(m*(br-lam))
assert R>1 and R<hi
print('PASS c1<0, c3>0 and strict residual-band margins at the tie')
for name,V in [('c1',c1),('c3',c3),('R',R),('upper',hi)]:print(name,M.nstr(M.mpf(V.mid),16))
