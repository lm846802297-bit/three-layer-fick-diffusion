from mpmath import iv
from fractions import Fraction
iv.dps=45
I=iv.mpf
pi=iv.pi
kappa=I(4); u=I(2);v=I(1)/5;z=I(4)/33
l1=(I(1)/2-z)*v; l2=z*u
def YL(k):
    return (iv.cos(k*l1)+kappa*v/k*iv.sin(k*l1))*iv.cos(k*l2) + (kappa*u/k*iv.cos(k*l1)-u/v*iv.sin(k*l1))*iv.sin(k*l2)
def r_eq(k,n):
    return k-(2*n-1)*pi-2*iv.atan2(kappa,k)
def root(n):
    lo,hi=Fraction(2*n-1)*Fraction('3.14159'),Fraction(2*n)*Fraction('3.14160')
    def val(x):return r_eq(I(x.numerator)/x.denominator,n)
    assert val(lo)<0 and val(hi)>0
    for _ in range(100):
        mid=(lo+hi)/2
        if val(mid)<0:lo=mid
        elif val(mid)>0:hi=mid
        else:raise ValueError('ambiguous bisection')
    return I([(I(lo.numerator)/lo.denominator).a,(I(hi.numerator)/hi.denominator).b])
def phase(k):
    d=k*l1-iv.atan2(kappa*v,k)
    j=0 if -pi/2<d<pi/2 else 1
    if j: assert pi/2<d<3*pi/2
    d-=j*pi
    return j*pi+iv.atan2(u/v*iv.sin(d),iv.cos(d))+k*l2
def left_eig(lo,hi,n):
    target=(I(n)-I(1)/2)*pi
    assert phase(I(lo))<target<phase(I(hi))
    return I([lo,hi])**2
assert (Fraction(1,2)-Fraction(4,33))*Fraction(1,25)+Fraction(4,33)*4==Fraction(1,2)
a=left_eig('3.77267640109','3.77267640111',1)
L2=left_eig('13.25724841026','13.25724841029',2)
N=43
ks=[root(n) for n in range(1,N+1)]
for n,k in enumerate(ks[:3],1):
    margin=phase(k)-(I(n)-I(1)/2)*pi
    assert (0<margin<pi) if n==1 else (-pi<margin<0)
assert 7*(I(1)/2-l1-l2)>1
b=ks[0]**2
m=ks[1]**2
q=1-m/L2

# Late endpoint: verify the exact-moment bound again using only these brackets.
def YR(k):return iv.cos(k/2)+kappa/k*iv.sin(k/2)
def dYR(k):return iv.sin(k/2)*(1+(kappa*kappa+2*kappa)/k**2)/(4*k)
g=I(183)
def ratio(eta):
    k=iv.sqrt(eta)
    return YL(k)/YR(k)*(1-eta/b)*(1-eta/m)/((1-eta/a)*(1-eta/L2))
Mg=ratio(g)
Mm=YL(ks[1])*(1-m/b)/(m*dYR(ks[1])*(1-m/a)*(1-m/L2))
late=I('0.0171')
assert (a+q*(m-a))*Mg*iv.exp(-(g-m)*late)<q*(m-a)*Mm
print('late bound',iv.ln((1+a/(q*(m-a)))*Mg/Mm)/(g-m))

# Early endpoint: on Re sqrt(s)>=delta/(2t0), |V(s)|<=25 exp(-delta Re sqrt(s)).
t0=I('0.0008'); delta=I(1)/2-l1-l2
Q=delta/(2*t0)
yb=((1-kappa/Q)-(1+kappa/Q)*iv.exp(-Q))/2
assert (1+u/v+kappa*(u+v)/Q)/yb<25
assert (a/b)*(1+b/Q**2)/(1-a/Q**2)<1
assert delta**2/4>I(3)/2*t0
B=25/pi*iv.exp(-delta**2/(4*t0))*(Q*iv.sqrt(pi/t0)+1/t0)
assert B<I("2.173")
assert B<a*(1-t0*B)
print('early upper f',B)

assert 11+I("8.8")/(3*pi)<12
assert 48*iv.sqrt(1+(4/(3*pi))**2)<60
pairs=[(k*k,a*(1/b-1/(k*k))*YL(k)/dYR(k)) for k in ks[1:]]
def bounds(lo,hi,order=0):
    ti=I([lo,hi])
    summ=sum(B*(-R)**order*iv.exp(-R*ti) for R,B in pairs)
    n=N+1; t=I(lo); exponent=1+2*order
    first=60*(2*pi*n)**exponent*iv.exp(-((2*n-1)*pi)**2*t)
    ratio=(I(n+1)/n)**exponent*iv.exp(-8*pi*pi*n*t)
    assert ratio<1
    tail=(first/(1-ratio)).b
    return summ+I([-tail,tail])
def exact_decimal(i):
    return f"0.{i:08d}"

def certify(lo,hi,sign,order):
    todo=[(lo,hi)]; count=0
    while todo:
        l,r=todo.pop()
        if sign*bounds(exact_decimal(l),exact_decimal(r),order)>0:
            count+=1
        else:
            assert r-l>1,(l,r,bounds(exact_decimal(l),exact_decimal(r),order))
            c=(l+r)//2
            todo.extend([(l,c),(c,r)])
    print('certified',lo,hi,'/1e8', 'sign',sign,'derivative order',order,'subintervals',count)
certify(80000,140000,-1,0)
certify(140000,170000,1,1)
certify(170000,1710000,1,0)
left,right=140000,170000
while right-left>1:
    mid=(left+right)//2
    x=exact_decimal(mid)
    H=bounds(x,x)
    if H<0:left=mid
    elif H>0:right=mid
    else:raise ValueError("could not bracket the root")
assert bounds(exact_decimal(left),exact_decimal(left))<0
assert bounds(exact_decimal(right),exact_decimal(right))>0
print('unique upward hazard threshold crossing in',exact_decimal(left),exact_decimal(right))

# Sign of the physical leading residue c1 using the full ground wave number.
def full_at(k):
    y=iv.cos(k*l1)+kappa*v/k*iv.sin(k*l1)
    q0=kappa*iv.cos(k*l1)-k/v*iv.sin(k*l1)
    left_y=YL(k)
    left_q=q0*iv.cos(k*l2)-(k/u)*y*iv.sin(k*l2)
    right_y=YR(k)
    right_q=kappa*iv.cos(k/2)-k*iv.sin(k/2)
    return left_y*right_q+right_y*left_q, left_y, right_y
left,right=I('2.41899445303'),I('2.41899445305')
assert full_at(left)[0]>0 and full_at(right)[0]<0
k=I(['2.41899445303','2.41899445305'])
assert k**2<a  # below the first split half eigenvalue; interlacing identifies lambda_1
_, yl, yr=full_at(k)
assert yl>0 and yr>yl
print('ground lambda_1',k**2)
print('ground right/left wall eigenfunction values',yl/yr)
print('c1>0: physical cumulative difference has exactly one strict reversal')
