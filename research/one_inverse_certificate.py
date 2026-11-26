import mpmath as m
m.mp.dps=65
m.iv.dps=55
I,M=m.iv,m.mp
def yR(q,A):return A.cos(q/2)+A.sin(q/2)/q
def lh(q,A):
    a=q*A.mpf('.455')/2
    y=A.cos(a)+A.sin(a)/(2*q)
    p=-2*A.sin(a)+A.cos(a)/q
    b=q*A.mpf('.09')
    return y*A.cos(b)+2*p*A.sin(b)
def phase(q,A):
    a=q*A.mpf('.455')/2-A.atan2(1,2*q)
    x,y=A.cos(a),A.sin(a)
    return a + A.atan2(3*x*y,x*x+4*y*y)+A.mpf('.09')*q
def isolate(fn,guess):
    k=M.findroot(lambda q:fn(q,M),guess)
    a,b=[M.nstr(k+d*M.mpf('1e-35'),60) for d in(-1,1)]
    assert fn(I.mpf(a),I)*fn(I.mpf(b),I)<0
    return I.mpf([a,b])
l=isolate(lh,(M.mpf('3.4'),M.mpf('3.6')))
assert abs(phase(l,I)-I.pi/2)<I.mpf('1e-32')
dL=I.mpf('.455')/2+I.mpf('.09')
assert 7*I.pi*(I.mpf('.5')-dL)>I.pi
for n in [1,2,3]:
    r=isolate(yR,((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    assert r>(2*n-1)*I.pi and r<2*n*I.pi
    g=phase(r,I)-(I.mpf(n)-I.mpf('.5'))*I.pi
    print('n',n,'phase gap at k_R',g)
    if n==1:
        assert l<r and g>0
        rho=(I.mpf('.455')/4+I.mpf('.045')*4+1)/I.mpf('1.5')
        bound=rho*(r/l)**2
        assert bound<I.mpf('.95')
        print('k_L1',l,'k_R1',r,'rho times p',bound)
    else:assert g<0
print('PASS exactly one reversed mode, all remaining ordered; no reversal, A_R(t)<.95 A_L(t)')
