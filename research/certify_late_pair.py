from certify_folds import calc,z,make
from upper import M,I,iv,mpf
import numpy as np

def s(x): return M.nstr(x,73)
N=120
def iso(K,q):
    rt=M.findroot(lambda k:calc(k,K,z(M,True),M),(mpf(float(q))-.001,mpf(float(q))+.001))
    lo,hi=[s(rt+d*mpf('1e-55')) for d in (-1,1)]
    assert calc(iv(lo),iv(s(K)),z(I,True),I)*calc(iv(hi),iv(s(K)),z(I,True),I)<0
    return lo,hi

def main():
    Z=z(I,True)
    q0=2.44
    def c1(K):
        r=M.findroot(lambda k:calc(k,K,z(M,True),M),q0)
        return calc(r,K,z(M,True),M,True)[0]
    k0=M.findroot(c1,(mpf('5.69'),mpf('5.71')))
    a,b=[mpf(s(k0+d*mpf('1e-48'))) for d in (-1,1)]
    K=iv([s(a),s(b)])
    for kap,sg in [(a,1),(b,-1)]:
        br=iso(kap,q0)
        w,_=calc(iv(br),iv(s(kap)),Z,I,True)
        assert sg*w>0
    print('PASS critical kappa in',s(a),s(b),flush=True)
    L=[];C=[]
    def floor(x):
        ll,rr=mpf(x.a),mpf(x.b)
        assert M.floor(ll)==M.floor(rr) and ll>M.floor(ll)
        return int(M.floor(ll))
    for q in np.sqrt(make(float(k0),z=float(z(M,True)),n=N)[1]):
        la,ra=iso(a,q); lb,rb=iso(b,q)
        k=iv([la,rb])
        w,ns=calc(k,K,Z,I,True)
        assert sum(floor(x)-floor(y) for x,y in ns)==len(L)
        assert not L or k**2>L[-1]
        L.append(k**2); C.append(w)
    C[0]=iv(0) # At a true critical parameter, exactly zero.
    print('PASS 120 ordered modes on critical enclosure',flush=True)
    # Same rational contour conditions as the previous benchmark,
    # uniformly on the critical kappa enclosure. (T=0.001)
    d1=(iv('.5')-Z)*iv('.04'); d2=Z*iv('1.1'); dr=iv('.5')
    dl=d1+d2; T=iv('.001'); qr=dr/(2*T); ql=dl/(2*T)
    gam=(1+25*iv('1.1'))*(1+1/iv('1.1'))/8
    cc=K/(2*gam); g=(dr**2-dl**2)/4
    U=[I.exp(-2*x*ql) for x in (d1,d2,dr)]
    ew=(1+U[0])*(1+U[1])*(1+U[2])-1
    assert ew<1 and (ew+U[2])/(1-ew)+4*I.sqrt(T)/(dl*I.sqrt(I.pi))+K*iv('.04')/ql<iv('.5')
    assert I.exp(-2*dr*qr)<iv(1)/3 and g>T/2
    assert 3*K*dl**2*I.sqrt(I.pi)/(cc*dr*I.sqrt(T))*I.exp(-g/T)<iv('.06')
    print('PASS 0<t<=.001: strict left lead, with uniform contour margins',flush=True)
    def tail(t,order):
        # kappa<6, lambda>=8*j*j, j>=120. For 0<=order<=2,
        # coefficient at minimal lambda bounded by 1000*j^3 and decreases
        # with lambda after multiplying by exp(-lambda*t). Ratio below 1/2.
        assert t>iv('.000998')
        return (200*N if order<=1 else 2000*N**3)*I.exp(-8*N*N*t)
    def point(t,d):
        v=-sum(c*(-l)**d*I.exp(-l*t) for l,c in zip(L,C))
        return v+iv([-1,1])*tail(t,d)
    P=12
    plan=[('.000999','.012',0,1),('.012','.013',1,-1),('.013','.5',0,-1)]
    for a,b,d,sg in plan:
        a,b=map(mpf,(a,b)); cells=0
        weights=[-c*(-l)**d for l,c in zip(L,C)]
        def check(a,b,depth=0):
            nonlocal cells
            mid=(iv(a)+iv(b))/2; rad=(iv(b)-iv(a))/2
            Q=[iv(0) for _ in range(P)]
            rem=iv(0)
            for l,w in zip(L,weights):
                f=w*I.exp(-l*mid)
                for j in range(P):
                    Q[j]+=f; f=f*(-l)/(j+1)
                rem+=abs(w)*I.exp(-l*iv(a))*(l*rad)**P/M.factorial(P)
            err=sum(abs(Q[j])*rad**j for j in range(1,P))+rem+tail(iv(a),d)
            if sg*Q[0]>err:cells+=1;return
            assert depth<26,(a,b,Q[0],err)
            m=(a+b)/2;check(a,m,depth+1);check(m,b,depth+1)
        check(a,b)
        print('PASS',s(a),s(b),d,sg,'slabs:',cells,flush=True)
    ta,tb=map(iv,('.01266184','.01266185'))
    assert point(ta,0)>0 and point(tb,0)<0
    print('root checks',point(ta,0),point(tb,0),flush=True)
    T=iv('.5')
    assert point(T,0)<0
    # Second derivative of exp(lambda1*t)*Delta; factor exp(-gap*t).
    lead=C[1]*(L[1]-L[0])**2
    hi=sum(abs(c)*(l-L[0])**2*I.exp(-(l-L[1])*T) for l,c in zip(L[2:],C[2:]))
    hi+=I.exp(L[1]*T)*tail(T,2) # (lambda_n-lambda_1)^2 < lambda_n^2
    assert hi<lead
    print('PASS tail second-derivative test at T=.5, lead',lead,'high-mode bound',hi,flush=True)
    assert C[1]>0
    # Also certifies Delta(T)<0 and a strict early slope sign on the root bracket.
if __name__=='__main__':main()
