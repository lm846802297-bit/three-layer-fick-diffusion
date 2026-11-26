"""Interval verification of the analytic all-time criterion at two examples."""
import mpmath as mp
M,I=mp.mp,mp.iv
M.dps=70; I.dps=60
def f(q,K,z,A,kind):
    c=A.mpf
    y,p=c(1),K/q
    left=((c('.5')-z,c('.02')),(z,c(2)))
    if kind=='right': ls=((c('.5'),c(1)),)
    else: ls=left+(((c('.5'),c(1)),) if kind=='full' else ())
    for h,r in ls:
        S,C=A.sin(q*h*r),A.cos(q*h*r)
        # y is positive throughout a segment of length < pi in phase
        # when both endpoint values are positive.
        y0=y
        y,p=y*C+r*p*S,p*C-y*S/r
        if kind=='check': yield y0,y,q*h*r
    if kind!='check': yield y if kind!='full' else p+K*y/q
def g(*a):return next(f(*a))
def bracket(K,z,kind,xa,xb):
    Q=M.findroot(lambda q:g(q,M.mpf(K),M.mpf(z),M,kind),(xa,xb),solver='anderson')
    edges=[M.nstr(Q+d*M.mpf('1e-34'),65) for d in [-1,1]]
    K,z=I.mpf(K),I.mpf(z)
    assert g(I.mpf(edges[0]),K,z,I,kind)*g(I.mpf(edges[1]),K,z,I,kind)<0
    return I.mpf(edges)
def check(Ks,zs):
    K,z=I.mpf(Ks),I.mpf(zs)
    h=I.mpf('.5')-z; v=I.mpf('.02')
    q=bracket(Ks,zs,'full',2,3)
    assert all(l>0 and r>0 and phase<I.pi for l,r,phase in f(q,K,z,I,'check'))
    # At the full eigenvalue the right segment is proportional (with
    # positive multiplier, matched at the midpoint) to the
    # right-wall solution, which has no node since q < pi.
    assert q<I.pi
    qa=bracket(Ks,zs,'left',3,5)
    qb=bracket(Ks,zs,'right',4,5)
    qm=bracket(Ks,zs,'right',9.5,11)
    a,b,m,L=qa**2,qb**2,qm**2,q**2
    assert qb>I.pi and qb<2*I.pi and qm>3*I.pi and qm<4*I.pi
    assert L<a and a<b
    t=qa*h*v-I.atan2(K*v,qa)
    c,s=I.cos(t),I.sin(t)
    phase=t+I.atan2(99*c*s,c*c+100*s*s)+2*z*qa
    assert phase>I.pi/4 and phase<3*I.pi/4 # first Dirichlet index
    for qr,n in [(qb,1),(qm,2)]:
        t=qr*h*v-I.atan2(K*v,qr)
        c,s=I.cos(t),I.sin(t)
        phase=t+I.atan2(99*c*s,c*c+100*s*s)+2*z*qr-(I.mpf(n)-I.mpf('.5'))*I.pi
        assert phase>0 if n==1 else phase<0
    assert 5*(I.mpf('.5')-h*v-2*z)>1
    rho=(h*v*v+4*z+1/K)/(I.mpf('.5')+1/K)
    R=g(q,K,z,I,'left')/g(q,K,z,I,'right')/rho
    lo=(rho*b-L)/(rho*(b-L)); hi=b*(m-L)/(rho*m*(b-L))
    assert rho<1 and lo<R and R<hi
    assert rho*b/a>1
    if Ks=='7' and zs=='.124': assert R<1
    print('PASS all-time left lead; kappa',Ks,'slow width',zs)
    for label,x in [('rho',rho),('lower',lo),('R',R),('upper',hi),('rho*p',rho*b/a)]:
        print(label,M.nstr(M.mpf(x.mid),14),'enclosure radius <',M.nstr(M.mpf(x.delta),3))
if __name__=='__main__':
    check('7','.124')
    check('7','.12')
