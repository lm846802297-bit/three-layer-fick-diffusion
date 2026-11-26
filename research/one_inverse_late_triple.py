import mpmath as mp
import numpy as np
from fold import make
M,I=mp.mp,mp.iv
M.dps=88;I.dps=82
def s(x):return M.nstr(x,85)
class AD:
 def __init__(self,x,d=None):self.x=I.mpf(x); self.d=[I.mpf(0)]*4 if d is None else d
 def __add__(a,b):
    if not isinstance(b,AD):b=AD(b)
    return AD(a.x+b.x,[x+y for x,y in zip(a.d,b.d)])
 __radd__=__add__
 def __neg__(a):return AD(-a.x,[-x for x in a.d])
 def __sub__(a,b):return a+-b
 def __rsub__(a,b):return -a+b
 def __mul__(a,b):
    if not isinstance(b,AD):b=AD(b)
    return AD(a.x*b.x,[x*b.x+a.x*y for x,y in zip(a.d,b.d)])
 __rmul__=__mul__
 def __truediv__(a,b):
    if not isinstance(b,AD):b=AD(b)
    return AD(a.x/b.x,[(x*b.x-a.x*y)/b.x**2 for x,y in zip(a.d,b.d)])
 def sqrt(a):
    r=I.sqrt(a.x);return AD(r,[d/(2*r) for d in a.d])
 def sin(a):return AD(I.sin(a.x),[I.cos(a.x)*d for d in a.d])
 def cos(a):return AD(I.cos(a.x),[-I.sin(a.x)*d for d in a.d])
def eq(args,ctx):
 h,K,l1,l2=args; P=ctx.mpf('0.5')
 if ctx is I and isinstance(h,AD):
    P=AD("0.5")
    root=lambda x:x.sqrt();sn=lambda x:x.sin();co=lambda x:x.cos()
 else:root,sn,co=ctx.sqrt,ctx.sin,ctx.cos
 def state(l,left):
    y,p=1,K; q=root(l)
    for width,res in ([(h,h/100),(P-h,P-h/100)] if left else [(P,P)]):
       L=root(width*res)*q; c=co(L); S=sn(L)/L
       y,p=c*y+res*S*p,-l*width*S*y+c*p
    return y,p
 a,b=state(l1,True);c,d=state(l1,False)
 return [a*d+b*c,a-c,state(l2,True)[0],state(l2,False)[0]]
def enclose(eps):
 v=M.findroot(lambda *a:eq(a,M),('0.188276062556','6.115183531695','6.668186310160','24.273455580416'))
 arr=[M.mpf(s(a)) for a in v]
 X=[I.mpf([s(a-eps),s(a+eps)]) for a in arr]
 mid=list(map(lambda a:I.mpf(s(a)),arr))
 Y=M.inverse(M.matrix(M.jacobian(lambda *a:eq(a,M),arr)))
 Y=[[I.mpf(s(Y[i,j])) for j in range(4)] for i in range(4)]
 A=[]
 for i,x in enumerate(X):
    row=[I.mpf(0)]*4;row[i]=I.mpf(1);A.append(AD(x,row))
 J=[a.d for a in eq(A,I)]
 F=eq(mid,I)
 for i in range(4):
    R=[int(i==j)-sum(Y[i][k]*J[k][j] for k in range(4)) for j in range(4)]
    V=-sum(Y[i][j]*F[j] for j in range(4))
    assert sum(abs(r) for r in R)<I.mpf('.001')
    assert abs(V+sum(R[j]*(X[j]-mid[j]) for j in range(4)))<I.mpf(s(eps))/2
 print('PASS Krawczyk and nonsingular Jacobian, radius',s(eps),flush=True)
 return arr,X
def calc(k,H,K,ctx,full=False):
 V=ctx.mpf
 p=K/k;y=V(1);n=V(0);ns=[]
 for i,(w,u) in enumerate([(H,V('.1')),(V('.5')-H,ctx.sqrt((V('.5')-H/100)/(V('.5')-H))),(V('.5'),V(1)),]):
    a,b,x=y,u*p,k*w*u
    co,si=ctx.cos(x),ctx.sin(x)
    if full:
       n+=(a*a+b*b)*w/2+(a*a-b*b)*ctx.sin(2*x)/(4*u*k)+a*b*si**2/(u*k)
       A=ctx.atan2(b,a);ns.append(((x-A-ctx.pi/2)/ctx.pi,(-A-ctx.pi/2)/ctx.pi))
    y,p=a*co+b*si,(b*co-a*si)/u
    if i==1:at=y
 return p+K*y/k if not full else (K*at*(1-y)/n/k**2,ns)
def main():
 arr,X=enclose(M.mpf('1e-58')); H,K=X[:2]
 h,k=arr[:2]; z=M.mpf('.5')-h;u=M.sqrt((M.mpf('.5')-h/100)/z)
 P=I.mpf
 q2=I.sqrt(X[3]); phase=q2*H/10
 assert phase<I.pi and I.cos(phase)+K/(10*q2)*I.sin(phase)>0
 assert q2*I.sqrt((P(".5")-H)*(P(".5")-H/100))<I.pi and q2/2<I.pi
 assert 1/((P('.5')-H/100)/(P('.5')-H))>6/I.pi**2 and K<8
 ell=H/10+I.sqrt((P('.5')-H)*(P('.5')-H/100))
 assert 13*(P('.5')-ell)>1
 def rf(q,Q,k):return Q.cos(q/2)+k*Q.sin(q/2)/q
 invs=[]
 for n in range(2,7):
    q=M.findroot(lambda k:rf(k,M,arr[1]),((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    lo,hi=[s(q+j*M.mpf('1e-53')) for j in [-1,1]]
    assert rf(P(lo),I,K)*rf(P(hi),I,K)<0
    q=P([lo,hi]); assert q>(2*n-1)*I.pi and q<2*n*I.pi
    u2=I.sqrt((P('.5')-H/100)/(P('.5')-H))
    a=q*H/10-I.atan2(K/10,q); x,y=I.cos(a),I.sin(a);r=10*u2
    d=a+I.atan2((r-1)*x*y,x*x+r*y*y)+q*(P('.5')-H)*u2-(P(n)-P('.5'))*I.pi
    if d>0:invs.append(n)
    else: assert d<0
 assert invs==[2]
 print('PASS critical half spectra: n=1 exact tie, n=2 inversion, all n>=3 normal',flush=True)
 L=[];C=[];N=140
 def fl(x):
    a,b=M.mpf(x.a),M.mpf(x.b)
    assert M.floor(a)==M.floor(b) and a>M.floor(a);return int(M.floor(a))
 for guess in np.sqrt(make(float(k),z=float(z),u=float(u),v=.1,n=N)[1]):
    q=M.findroot(lambda x:calc(x,h,k,M),(M.mpf(float(guess))-.001,M.mpf(float(guess))+.001))
    a,b=[s(q+j*M.mpf('1e-53')) for j in [-1,1]]
    assert calc(P(a),H,K,I)*calc(P(b),H,K,I)<0
    Q=P([a,b]); Xq=Q**2
    c,ns=calc(Q,H,K,I,True)
    # The second mode has an exact node at layer boundary x=.5; its
    # index is checked analytically from the two positive half-solutions.
    if len(L)!=1:assert sum(fl(t)-fl(v) for t,v in ns)==len(L)
    else:assert Xq<X[3]+P('1e-51') and Xq>X[3]-P('1e-51')
    assert not L or Xq>L[-1]
    L.append(Xq); C.append(c)
 assert len(L)==N
 assert C[2]<0;C[:2]=[P(0)]*2
 print('PASS modal isolation; original-convention c3',C[2],flush=True)
 d1=H/10;d2=I.sqrt((P('.5')-H)*(P('.5')-H/100));dr=P('.5');dl=d1+d2;T=P('.001')
 Q=dl/(2*T); qr=dr/(2*T)
 Z2=(P('.5')-H)/d2
 gm=(1+10/Z2)*(1+Z2)/8; g=(dr**2-dl**2)/4; coef=K/(2*gm)
 U=[I.exp(-2*d*Q) for d in (d1,d2,dr)]
 E=(1+U[0])*(1+U[1])*(1+U[2])-1
 assert E<1 and (E+U[2])/(1-E)+4*I.sqrt(T)/(dl*I.sqrt(I.pi))+K/(10*Q)<P('.5')
 assert I.exp(-2*dr*qr)<P(1)/3 and g>T/2
 assert 3*K*dl**2*I.sqrt(I.pi)/(coef*dr*I.sqrt(T))*I.exp(-g/T)<1
 print('PASS strict left lead 0<t<=.001',flush=True)
 # D_min>.625 -> lambda_n>6(n-1)^2. kappa<8;
 # at lambda>=6*j^2, j>=140, flux bound 16*(1+2.53*sqrt(lambda))
 # gives <150*j at equality, decreases after multiplication by exp(-lambda*t).
 # Adjacent tail bounds have ratio <1/2 for t>=.000999.
 def tail(t):return 300*N*I.exp(-6*N*N*t)
 def check(a,b,d=0):
    nonlocal cells
    mid=(P(a)+P(b))/2;rad=(P(b)-P(a))/2; Q=[P(0)]*12; rem=P(0)
    for l,c in zip(L[2:],C[2:]):
       f=-c*I.exp(-l*mid)
       for j in range(12):
           Q[j]+=f;f=f*(-l)/(j+1)
       rem+=abs(c)*I.exp(-l*P(a))*(l*rad)**12/M.factorial(12)
    B=sum(abs(Q[j])*rad**j for j in range(1,12))+rem+tail(P(a))
    if Q[0]>B:cells+=1;return
    assert d<26
    m=(a+b)/2;check(a,m,d+1);check(m,b,d+1)
 cells=0;check(M.mpf('.000999'),M.mpf('.3'))
 print('PASS positivity .000999 through .3, slabs',cells,flush=True)
 t=P('.3')
 hi=sum(abs(c)*I.exp(-(l-L[2])*t) for l,c in zip(L[3:],C[3:]))+I.exp(L[2]*t)*tail(t)
 assert hi<-C[2]
 print('PASS positivity beyond .3 by third mode dominance, tail:',hi,flush=True)
 # Apply three weighted derivatives to exp(lambda3*t)*Delta;
 # they remove d, c1, c2. The upper bound is decreasing for t>=T.
 def wt(l):return l*(l-L[0])*(l-L[1])
 lead=-C[2]*wt(L[2])
 hi=sum(abs(c)*wt(l)*I.exp(-(l-L[2])*t) for l,c in zip(L[3:],C[3:]))
 # For omitted terms use wt(lambda)<lambda^3; pointwise coefficient
 # 16*(1+2.53*sqrt(lambda))*lambda^2. At lambda=6*j^2 this
 # is <4000*j^5, and exponential monotonicity and ratio give the bound.
 hi+=I.exp(L[2]*t)*8000*N**5*I.exp(-6*N*N*t)
 assert hi<lead
 print('PASS tail weighted-three-derivative bound, lead',lead,'remaining',hi,flush=True)
if __name__=='__main__':main()
