"""Independent interval check at kappa=8 for the static tie example."""
from fractions import Fraction as Q
from mpmath import iv
I=iv.mpf; iv.dps=35; p=iv.pi; k=I(8); u=I(2); v=I(1)/5; z=I(4)/33
l=(I(1)/2-z)*v; m=z*u
assert Q(25,66)/25+Q(4,33)*4==Q(1,2)
rs=[("5.14086312067","5.14086312068"),("10.70806368234","10.70806368235"),("16.60585836519","16.60585836520")]
def phase(q):
    x=q*l-iv.atan2(k*v,q)
    if -p/2<x<p/2:return iv.atan2(u/v*iv.sin(x),iv.cos(x))+q*m
    assert p/2<x<3*p/2
    x-=p
    return p+iv.atan2(u/v*iv.sin(x),iv.cos(x))+q*m
for j,(a,b) in enumerate(rs,1):
    eq=lambda q:q-(2*j-1)*p-2*iv.atan2(k,q)
    assert eq(I(a))<0<eq(I(b))
    d=phase(I([a,b]))-(I(j)-I(1)/2)*p
    assert (0<d<p) if j==1 else (-p<d<0)
    print("half phase margin",j,d)
assert 7*(I(1)/2-l-m)>1
def W(q):
    x,y,w=q*l,q*m,q/2
    L=iv.cos(x)+k*v/q*iv.sin(x); A=k*iv.cos(x)-q/v*iv.sin(x)
    Y=L*iv.cos(y)+u/q*A*iv.sin(y); QL=A*iv.cos(y)-q/u*L*iv.sin(y)
    R=iv.cos(w)+k/q*iv.sin(w); S=k*iv.cos(w)-q*iv.sin(w)
    return Y*S+QL*R,Y,R
a,b="2.96867131114","2.96867131115"
assert W(I(a))[0]>0>W(I(b))[0]
assert phase(I([a,b]))<p/2
assert I(b)<I(rs[0][0])
_,Y,R=W(I([a,b])); assert Y>R>0
print("ground amplitude ratio right/left",Y/R)
print("c1 < 0; by the single-zero theorem, there is no capture reversal")
