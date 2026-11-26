"""Verify a genuine first-only reversal point with d=c1=0 and full-rank tuning."""
from late_triple import AD
import mpmath as m
M,I=m.mp,m.iv
M.dps=75;I.dps=65
P=I.mpf
def s(x):return M.nstr(x,72)
def eq(q,K,A,where=False):
    ad=isinstance(q,AD)
    Q=AD if ad else A.mpf
    zm=Q('.5')*(1-Q('.02')*Q('.02'))/(4-Q('.02')*Q('.02'))
    u=Q(2);v=Q('.02')
    y,p=Q(1),K/q
    states=[]
    for w,r in [(Q('.5')-zm,v),(zm,u),(Q('.5'),Q(1))]:
        L=q*w*r
        C,S=(L.cos(),L.sin()) if ad else (A.cos(L),A.sin(L))
        y,p=y*C+r*p*S,-y*S/r+p*C
        states.append((y,L))
    return states if where else (p+K*y/q,y-1)
cen=M.findroot(lambda q,k:eq(q,k,M),('2.8933','7.07884'))
eps=M.mpf('1e-12')
box=[P([s(x-eps),s(x+eps)]) for x in cen];mid=[P(s(x)) for x in cen]
Y0=M.inverse(M.matrix(M.jacobian(lambda q,k:eq(q,k,M),cen)))
Y=[[P(s(Y0[i,j])) for j in range(2)] for i in range(2)]
V=eq(*mid,I)
q=AD(box[0],[P(1),P(0),P(0),P(0)])
k=AD(box[1],[P(0),P(1),P(0),P(0)])
J=[h.d for h in eq(q,k,I)]
for i in range(2):
    R=[int(i==j)-sum(Y[i][l]*J[l][j] for l in range(2)) for j in range(2)]
    assert sum(abs(r) for r in R)<P('.001')
    assert abs(-sum(Y[i][j]*V[j] for j in range(2))+sum(R[j]*(box[j]-mid[j]) for j in range(2)))<P(s(eps))/2
print('PASS contraction existence and nonzero Jacobian',s(cen[0]),s(cen[1]))
q,K=box
assert all(y>0 and x<I.pi for y,x in eq(q,K,I,True))
Z=P('.5')*(1-P('.02')**2)/(4-P('.02')**2)
h=P('.5')-Z
gap=P('.5')-h*P('.02')-2*Z
assert 5*gap>1
def yr(q,k,A):return A.cos(q/2)+k*A.sin(q/2)/q
deltas=[]
for n in [1,2]:
    rr=M.findroot(lambda k:yr(k,cen[1],M),((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    a,b=[P(s(rr+d*eps*10)) for d in (-1,1)]
    assert a>(2*n-1)*I.pi and b<2*n*I.pi
    assert yr(a,K,I)*yr(b,K,I)<0
    Q=P([M.mpf(a.a),M.mpf(b.b)])
    A=Q*h*P('.02')-I.atan2(K*P('.02'),Q)
    c,t=I.cos(A),I.sin(A)
    d=A+I.atan2(99*c*t,c*c+100*t*t)+2*Z*Q-(n-.5)*I.pi
    deltas.append(d)
assert deltas[0]>0 and deltas[1]<0
print('PASS ground state index n=1 and unique half-spectrum inversion n=1',deltas)
