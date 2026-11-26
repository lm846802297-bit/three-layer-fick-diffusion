"""Certify the weighted derivative at the double-cancellation geometry."""
import first_inversion_weighted_late_pair as base
M,I,P,s=base.M,base.I,base.P,base.s
z,k,a,b=base.box
N=145
L=[]; C=[]
roots=base.mode_factory(float(base.us),float(base.vs),float(base.X0[0]),False)[0]
def res(q,z,k,A):
    y,p=base.state(q,z,k,A,'full')
    return p+k*y/q
for i in range(N):
    if i==1:q=a
    elif i==2:q=b
    else:
        x=roots(i,float(base.X0[1]))
        p=M.findroot(lambda q:res(q,*base.X0[:2],M),(x-.001,x+.001))
        A,B=[P(s(p+d*M.mpf('1e-48'))) for d in (-1,1)]
        assert res(A,z,k,I)*res(B,z,k,I)<0
        q=P([M.mpf(A.a),M.mpf(B.b)])
    y,g,x,J,c=base.state(q,z,k,I,'modal')
    if i==1 or i==2:c=P(0) # exact analytic cancellation, not floating subtraction
    else:assert base.index(J)==i
    assert not L or q*q>L[-1]
    L.append(q*q); C.append(c)
ell=P('.5')+(P('.5')-z)*P(base.vs)+z*P(base.us)
assert ell<1 and k<4 and P(base.us)<5
# Full phase < q*ell + pi, so lambda_n > 9(n-2)^2. For a normalized
# eigenfunction sup|phi|^2 <= 1+2 sqrt(lambda/Dmin), hence |lambda c| < 8(1+10sqrt(lambda)).
def tail(t):
    j=N-1; h=9*j*j
    assert h*t>P('1.5')
    assert (P(j+1)/j)**3*I.exp(-9*(2*j+1)*t)<P('.5')
    return 16*(1+10*I.sqrt(P(h)))*h*I.exp(-(h-L[3])*t)
G=[P(0)]+[l-L[3] for l in L[4:]]
W=[c*l*(L[0]-l) for l,c in zip(L[3:],C[3:])]
assert W[0]>0 and all(g>0 for g in G[1:])
assert W[0]>sum(abs(w)*I.exp(-g*P('.10')) for w,g in zip(W[1:],G[1:]))+tail(P('.10'))
print('PASS e^(lambda4-lambda2)t W is positive for t>=.1')
def certify(aa,bb):
    stack=[(M.mpf(aa),M.mpf(bb),0)]; count=0
    while stack:
        a,b,d=stack.pop();m=(a+b)/2;rad=(P(b)-P(a))/2
        M0=12; A=[P(0)]*M0;Rem=P(0)
        for g,w in zip(G,W):
            e=w*I.exp(-g*P(m))
            for j in range(M0):
                A[j]+=e; e*=(-g)/(j+1)
            Rem+=abs(w)*I.exp(-g*P(a))*(g*rad)**M0/M.factorial(M0)
        if A[0]>Rem+sum(abs(A[j])*rad**j for j in range(1,M0))+tail(P(a)):
            count+=1
        else:
            assert d<21,(a,b,A[0])
            stack.extend([(a,m,d+1),(m,b,d+1)])
    print('PASS complete slab',aa,bb,'pieces',count,flush=True)
for a,b in [('.002','.003'),('.003','.005'),('.005','.010'),('.010','.10')]:certify(a,b)
# Apply one more deflation to W': constant and mode 3 drop out exactly.
t=P('.10')
alpha=L[2]-L[1]
g4=L[3]-L[1]
lead=g4*(g4-alpha)*W[0]
err=sum((l-L[1])*(l-L[2])*abs(w)*I.exp(-g*t) for l,g,w in zip(L[4:],G[1:],W[1:]))
j=N-1; h=9*j*j
assert h*t>P('3.5')
assert (P(j+1)/j)**7*I.exp(-9*(2*j+1)*t)<P('.5')
err+=16*(1+10*I.sqrt(P(h)))*h**3*I.exp(-(h-L[3])*t)
assert lead>err
print('PASS strict margin for two-mode-deflated derivative on t>=.1, ratio',M.nstr(M.mpf((err/lead).mid),13))
