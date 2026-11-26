"""Exact rational certificate of >=3 lead reversals in a layered Robin interval."""
from fractions import Fraction as F
from math import factorial
from decimal import Decimal, localcontext

def point(x):
    return (F(x), F(x))
def add(*terms):
    return tuple(sum(y,F(0)) for y in zip(*terms))
def mul(x,y):
    # Used only for nonnegative quantities.
    assert x[0]>=0 and y[0]>=0
    return (x[0]*y[0],x[1]*y[1])
def h(x, parity, N=80):
    # Positive Taylor series; the omitted tail is bounded geometrically.
    assert x>=0 and parity in (0,1)
    l=sum(x**j/factorial(j) for j in range(parity,2*N+parity,2))
    j=2*N+parity
    r=x*x/((j+1)*(j+2))
    assert r<1
    return (l,l+(x**j/factorial(j))/(1-r))
def layer(y,p,d,u,q):
    c,s=h(d*q,0),h(d*q,1)
    return add(mul(c,y),mul(mul(point(u),s),p)),add(mul(mul(point(1/u),s),y),mul(c,p))
def fmt(x):
    with localcontext() as ctx:
        ctx.prec=18
        return str(Decimal(x.numerator)/Decimal(x.denominator))
def main():
    u,v,z,k=F('1.1'),F('.04'),F('.4132'),F('6.5')
    w=F(1,2)-z
    print('Dfast, Dslow, left optical length:',1/v**2,1/u**2,w*v+z*u)
    RL=w*v*v+z*u*u+1/k; RR=F(1,2)+1/k
    steady=(RR-RL)/(RR+RL)
    assert steady<0 and w*v+z*u<F(1,2)
    print('Delta infinity:',steady,fmt(steady))
    # Cross-check exponential-transform signs. N is right minus left shape,
    # not a probability; the denominator of the transform is positive.
    for q, low, high in [(F(1,10),'-.0008','-.0006'),(F(1),'0.0008','0.0009'),
                         (F(10),'-4','-1'),(F(25),'1000','20000')]:
        y0,p0=point(1),point(k/q)
        left,_=layer(*layer(y0,p0,w*v,v,q),z*u,u,q)
        right,_=layer(y0,p0,F(1,2),F(1),q)
        lo,hi=right[0]-left[1],right[1]-left[0]
        print('s =',q*q,'; N in',fmt(lo),fmt(hi))
        assert F(low)<lo<=hi<F(high)
    print('PASS: four strictly alternating signs in the real Laplace transform.')
if __name__=='__main__':
    main()
