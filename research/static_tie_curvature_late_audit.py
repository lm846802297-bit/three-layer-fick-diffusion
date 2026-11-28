from mpmath import iv
iv.dps=50
I=iv.mpf
pi=iv.pi
a=I(['3.77267640109','3.77267640111'])**2
l=I(5)/66
n=I(8)/33
def root(j):
    def fun(x):return x-(2*j-1)*pi-2*iv.atan2(4,x)
    lo=I('3.14159265358979323846264338327950288')*(2*j-1)
    hi=I('3.14159265358979323846264338327950289')*2*j
    assert fun(lo)<0<fun(hi)
    for _ in range(96):
        mid=(lo+hi)/2
        y=fun(mid)
        if y.b<0:lo=mid
        elif y.a>0:hi=mid
        else:raise ArithmeticError
    return I([lo.a,hi.b])
b=root(1)**2
N=18
D=[]
for j in range(2,N+1):
    k=root(j); R=k*k
    YL=iv.cos(n*k)*iv.cos(l*k)-10*iv.sin(n*k)*iv.sin(l*k)+4/k*(iv.cos(n*k)*iv.sin(l*k)/5+2*iv.sin(n*k)*iv.cos(l*k))
    dy=iv.sin(k/2)*(1+24/R)/(4*k)
    coef=(a/b)*(R-b)/(R-a)*YL/dy
    D.append((R,coef))
base=D[0][0]+D[1][0]
pos=D[0][1]*D[1][1]*(D[1][0]-D[0][0])**2
assert pos>0
assert 48*(1+I(16)/(3*pi)**2)**I('.5')<60
assert I(11)+I('8.8')/(3*pi)<12
def check(t):
    t=I(t)
    bad=I(0)
    for i,(L,A) in enumerate(D):
        for j,(R,B) in enumerate(D[i+1:],i+1):
            if i==0 and j==1:continue
            assert L+R>base
            q=A*B*(R-L)**2
            if not q.a>=0:bad+=(-q.a)*iv.exp(-(L+R-base)*t)
    def bound(j,p):return 60*(2*j*pi)**p*iv.exp(-((2*j-1)*pi)**2*t)
    def tail(j,p):
        ratio=(I(j+1)/j)**p*iv.exp(-8*j*pi*pi*t)
        assert ratio<1
        return bound(j,p)/(1-ratio)
    abs1=sum(abs(A)*iv.exp(-R*t) for R,A in D)+tail(N+1,1)
    error=iv.exp(base*t)*abs1*tail(N+1,5)
    return (pos-bad-error)/pos
m=check('.012')
assert m>0
print('positive leading coefficient',pos)
print('normalized margin at and beyond 0.012',m)
