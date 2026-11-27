"""Interval check for the uniform-right late hazard bound. Run with mpmath."""
from fractions import Fraction as Q
from mpmath import iv

iv.dps = 42
I = iv.mpf
v, u, z, kap = I(1)/5, I(2), I(4)/33, I(4)
pi = iv.pi
ell = (I(1)/2-z)*v + z*u
assert (Q(1,2)-Q(4,33))*Q(1,5)**2+Q(4,33)*2**2 == Q(1,2)
assert 7*(I(1)/2-ell)>1  # n >= 4: existing optical phase bound

def left_phase(k):
    theta = k*(I(1)/2-z)*v-iv.atan2(kap*v,k)
    # Only the first two possible branches are needed at the enclosed points.
    if theta > -pi/2 and theta < pi/2:
        j=0
    elif theta > pi/2 and theta < 3*pi/2:
        j=1
    else:
        raise ValueError("branch not certified")
    x=theta-j*pi
    return j*pi+iv.atan2((u/v)*iv.sin(x),iv.cos(x))+k*z*u

def right_eq(k,n):
    return k-(2*n-1)*pi-2*iv.atan2(kap,k)

def bracket(lo,hi,residual):
    assert residual(I(lo))<0 and residual(I(hi))>0
    return I([lo,hi])

rw=[bracket(lo,hi,lambda k,n=n:right_eq(k,n))
    for n,(lo,hi) in enumerate([
        ("4.57785945620","4.57785945622"),
        ("10.17397018820","10.17397018821"),
        ("16.19232720644","16.19232720646")],1)]
for n,k in enumerate(rw,1):
    margin=left_phase(k)-(I(n)-I(1)/2)*pi
    assert (0<margin<pi) if n==1 else (-pi<margin<0)
a=bracket("3.77267640109","3.77267640111",lambda k:left_phase(k)-pi/2)**2
L2=bracket("13.25724841026","13.25724841029",lambda k:left_phase(k)-3*pi/2)**2
b,m,next_rate=(k**2 for k in rw)
assert a<b and m<L2
g=I(183)
assert m<g<next_rate
q=1-m/L2
w=iv.sqrt(g)
Yr=(iv.cos(w/2)+(kap/w)*iv.sin(w/2))/(1+kap/2)
M=(1-g/b)*(1-g/m)/Yr
prefactor=M*(1+a/(q*(m-a)))
assert prefactor*iv.exp(-(g-m)*I("0.041"))<1
print("Only first half mode inverted; d=0 by exact rational geometry.")
for name,value in [("a",a),("R2",m),("L2",L2),("q",q),("suffix moment upper bound",M),
                   ("cutoff upper bound",iv.ln(prefactor)/(g-m))]:
    print(name,value)
print("h(t)>a certified for all t >= 0.041; no claim below cutoff.")
