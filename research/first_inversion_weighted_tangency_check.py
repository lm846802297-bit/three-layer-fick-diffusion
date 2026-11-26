"""Interval check of closest floating-point extremum and late derivative bounds."""
import mpmath as mp
from scan_critical import mode_factory
I,M=mp.iv,mp.mp
I.dps=64; M.dps=74
P=I.mpf

def run(us,vs,zs,Ks,N):
    def setup(q,A,cap,phase=False):
        f=A.mpf; u,v,z,K=map(f,(us,vs,zs,Ks)); y,p=f(1),K/q
        ns=[]; norm=f(0); angle=-A.atan2(K*v,q)
        for j,(h,r) in enumerate(((f('.5')-z,v),(z,u),(f('.5'),f(1)))):
            if j==2 and phase:return angle
            a,b,x=y,r*p,q*h*r; S,C=A.sin(x),A.cos(x); angle+=x
            if j==0:
                c,s=A.cos(angle),A.sin(angle)
                angle+=A.atan2((u/v-1)*c*s,c*c+(u/v)*s*s)
            if cap:
                norm+=h*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*q)
                ang=A.atan2(b,a); ns.append(((x-ang-A.pi/2)/A.pi,(-ang-A.pi/2)/A.pi))
            y,p=a*C+b*S,(b*C-a*S)/r
            if j==1:mid=y
        if cap:return K*mid*(1-y)/(norm*q*q),ns
        return p+K*y/q
    uu,vv,zz,kk=map(P,(us,vs,zs,Ks))
    gap=P('.5')-(P('.5')-zz)*vv-zz*uu
    assert gap>0 and (P('.5')-zz)*vv**2+zz*uu**2<P('.5')
    def root(g,q):
        center=M.findroot(g,(q-M.mpf('.00001'),q+M.mpf('.00001')))
        return [M.nstr(center+d*M.mpf('1e-42'),70) for d in (-1,1)]
    for n in range(1,110):
        if n>1 and (2*n-1)*gap>=1:break
        h=lambda x,A:x/2-A.atan2(A.mpf(Ks),x)-(n-A.mpf('.5'))*A.pi
        q=M.findroot(lambda x:h(x,M),((2*n-1)*M.pi,(2*n)*M.pi))
        a,b=root(lambda x:h(x,M),q)
        assert h(P(a),I)*h(P(b),I)<0
        assert (2*n-1)*I.pi<P(a) and P(b)<2*n*I.pi
        v=setup(P([a,b]),I,False,True)-(n-P('.5'))*I.pi
        assert v>0 if n==1 else v<0
    else: raise AssertionError('gap threshold not reached')
    q0=mode_factory(float(us),float(vs),float(zs),False)[0]
    L=[]; C=[]
    def idx(q):
        l,r=M.mpf(q.a),M.mpf(q.b)
        assert M.floor(l)==M.floor(r) and l>M.floor(l)
        return int(M.floor(l))
    for n in range(N):
        a,b=root(lambda q:setup(q,M,False),M.mpf(q0(n,float(Ks))))
        assert setup(P(a),I,False)*setup(P(b),I,False)<0
        Q=P([a,b]); c,ns=setup(Q,I,True)
        assert sum(idx(x)-idx(y) for x,y in ns)==n
        L.append(Q*Q); C.append(c)
    W=[c*l*(L[0]-l) for l,c in zip(L[1:],C[1:])]
    assert W[0]>0 and kk<6 and uu<=2
    # ell_total<1, with two interface phase jumps each <pi/2:
    # lambda_n>9(n-2)^2. |lambda*c_n|<=12(1+4*sqrt(lambda)).
    def tail(t,d,shift):
        j=N-1; k=9*j*j
        assert (k*t)>P('3.5')
        assert (P(j+1)/j)**7*I.exp(-9*(2*j+1)*t)<P('.5')
        return 24*(1+4*I.sqrt(P(k)))*k**(1+d)*I.exp(-(k-shift)*t)
    return L,W,tail

if __name__=='__main__':
    L,W,tail=run('1.1155035225034142','0.7823106852553852','0.30678997573023653','0.06667564120733163',115)
    def val(t,d):
        return sum(w*(L[1]-l)**d*I.exp((L[1]-l)*t) for w,l in zip(W,L[1:]))+P([-1,1])*tail(t,d,L[1])
    a,b='0.0064015615','0.0064015627'
    B=P([a,b])
    assert val(P(a),1)>0 and val(P(b),1)<0
    assert val(B,2)<0 and val(B,0)>P('.00160')
    print('PASS closest screened maximum is strictly above zero on bracket',a,b,val(B,0))
    for case in [('2','.02','.12','5.0902282130'),('1.9','.03','.116','5')]:
        L,W,tail=run(*case,120)
        t=P('.03')
        lead=-(L[2]-L[1])*W[1]
        err=sum((l-L[1])*abs(w)*I.exp((L[2]-l)*t) for l,w in zip(L[3:],W[2:]))+tail(t,1,L[2])
        assert lead>err
        print('PASS W strictly increasing for t>=.03:',case,'error/lead =',M.nstr(M.mpf((err/lead).mid),12))
