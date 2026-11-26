import mpmath as m
m.mp.dps=65;m.iv.dps=55
I=m.iv; M=m.mp
# Left: 0.49 of D=4, 0.01 of D=1/4. Right: 0.5 of D=1.
# Identical Robin coefficient kappa=1.
def yR(q,A):return A.cos(q/2)+A.sin(q/2)/q
dL=I.mpf('.49')/2+I.mpf('.01')*2
assert I.mpf('.49')/4+I.mpf('.01')*4<I.mpf('.5')
assert 5*I.pi*(I.mpf('.5')-dL)>I.pi
for n in [1,2]:
    k=M.findroot(lambda q:yR(q,M),((2*n-1)*M.pi+.01,2*n*M.pi-.01))
    a,b=[M.nstr(k+d*M.mpf('1e-35'),60) for d in (-1,1)]
    assert I.mpf(a)>(2*n-1)*I.pi and I.mpf(b)<2*n*I.pi
    assert yR(I.mpf(a),I)*yR(I.mpf(b),I)<0
    q=I.mpf([a,b])
    phase=q*I.mpf('.49')/2-I.atan2(1,2*q)
    x,y=I.cos(phase),I.sin(phase)
    rotated=I.atan2(3*x*y,x*x+4*y*y)
    gap=phase+rotated+I.mpf('.02')*q-(I.mpf(n)-I.mpf('.5'))*I.pi
    assert gap<0
    print('n=',n,'q_R=',a,b,'left phase minus required phase=',gap)
print('PASS all n>=3 from optical-length phase bound')
