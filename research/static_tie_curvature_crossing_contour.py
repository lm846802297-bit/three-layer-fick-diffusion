from mpmath import iv
from sys import argv
iv.dps=55
I=iv.mpf
pi=iv.pi
T0=I('.010'); T1=I('.012'); Q=int(argv[1]) if len(argv)>1 else 24; A=12; M=260
assert Q in (24,26)
lo,hi=Q-A,Q+A
print('Q =',Q,' strip =',A,' N =',M,flush=True)
k=I(4); d=I(2)/11; u=I(2); v=I(1)/5
l=I(25)/66*v; n=I(4)/33*u
def phase(x):
    w=x*l-iv.atan2(k*v,x)
    assert -pi/2<w<pi/2
    return iv.atan2(10*iv.sin(w),iv.cos(w))+x*n
x,y=map(I,['3.77267640109','3.77267640111'])
assert phase(x)<pi/2<phase(y)
for _ in range(95):
    m=(x+y)/2
    z=phase(m)-pi/2
    if z.b<0:x=m
    elif z.a>0:y=m
    else:raise ArithmeticError('wide left root')
a=I([x.a,y.b])**2
def re(x):return x-pi-2*iv.atan2(k,x)
x,y=I('4'),I('5')
assert re(x)<0<re(y)
for _ in range(110):
    m=(x+y)/2
    if re(m)<0:x=m
    elif re(m)>0:y=m
    else:raise ArithmeticError('wide root')
b=I([x.a,y.b])**2
assert 0<a<b and b-a<10
def chsh(z):
    w=iv.exp(z); W=1/w
    return (w+W)/2,(w-W)/2
def transform(q):
    c,s=chsh(l*q); C,S=chsh(n*q); U,V=chsh(q/2)
    return ((C*c+10*S*s+k/q*(v*C*s+u*S*c))/(U+k/q*V))*((1+q*q/b)/(1+q*q/a))
# Holomorphy and uniform transform majorant for Re q in [Q-A,Q+A].
assert 2*(11+I('8.8')/lo)/(1-4/I(lo)-(1+4/I(lo))*iv.exp(-lo))*(1+10/I(lo**2))<40
# Poisson alias error, integrating the strip boundary majorant.
# Integral_0^infty y^j exp(-T0*y^2)dy, by recursion.
from math import comb, factorial
mom=[iv.sqrt(pi/T0)/2,1/(2*T0)]
for j in range(2,8):mom.append((j-1)*mom[j-2]/(2*T0))
U=80/pi*iv.exp(T1*hi**2)*sum(comb(7,j)*hi**(7-j)*mom[j] for j in range(8))
r=iv.exp(-2*pi*A)
alias=2*U*r/(1-r)
# Omitted trapezoid terms, both signs of y.
s=M+1
r=(I(s+1)/s)**7*iv.exp(-T0*(2*s+1))
assert r<1
tail=80/pi*iv.exp(T1*Q**2)*(Q+s)**7*iv.exp(-T0*s*s)/(1-r)
E=(alias+tail).b
print('quadrature error',E,flush=True)
assert E<I('1e-10')
err=I([-E,E])
st=I('.00001'); w=st/2; degree=4
terms=[]; rb=[I(0)]*4
for j in range(M+1):
    q=iv.mpc(Q,j); p=q*q
    z=q/pi*transform(q)*(2 if j else 1)
    amp=iv.exp(p.real*(T1 if p.real>=0 else T0))
    for h in range(4):
        rb[h]+=abs(z)*amp*abs(p)**(h+degree+1)*w**(degree+1)/factorial(degree+1)
    terms.append([z*iv.exp(p*(T0+w)),iv.exp(p*st),[p**z for z in range(degree+4)]])
# Use separate powers of the symmetric cell interval to preserve even-power positivity.
Wp=[I(1)]+[I([0,(w**k).b]) if k%2==0 else I([-(w**k).b,(w**k).b]) for k in range(1,degree+1)]
def bounds(vals,wp,rem):
    return [sum(vals[h+k]*wp[k]/factorial(k) for k in range(len(wp)))+rem[h]+err for h in range(4)]
def metrics(jets):
    f,g,j,k=jets
    assert f>0
    v=g/f
    return j/f-v**2,k/f-3*v*j/f+2*v**3
rem=[I([-x.b,x.b]) for x in rb]
minimum=I('1e20')
for t in range(200):
    vals=[I(0)]*(degree+4)
    for row in terms:
        z,ratio,powers=row
        for m,p in enumerate(powers):vals[m]+=(z*p).real
        row[0]=z*ratio
    f=bounds(vals,Wp,rem)
    c,dc=metrics(f)
    assert dc>0,(t,dc)
    minimum=min(minimum,dc.a)
print('continuous 200 cells; minimum curvature derivative',minimum)
def at(t):
    t=I(t)
    vals=[I(0)]*4
    for j in range(M+1):
        q=iv.mpc(Q,j); p=q*q
        z=q/pi*transform(q)*(2 if j else 1)*iv.exp(p*t)
        for h in range(4): vals[h]+=(z*p**h).real
    return metrics([x+err for x in vals])[0]
from decimal import Decimal
left=Decimal('.010'); right=Decimal('.012')
assert at(str(left))<0<at(str(right))
print('curvature endpoints',at(str(left)),at(str(right)))
while right-left>Decimal('1e-10'):
    for ratio in (Decimal('.5'),Decimal('.25'),Decimal('.75')):
        t=left+(right-left)*ratio
        c=at(str(t))
        if c.b<0:
            left=t
            break
        if c.a>0:
            right=t
            break
    else:
        raise ArithmeticError('root interval too wide')
print('unique positive-going root:',left,right)
