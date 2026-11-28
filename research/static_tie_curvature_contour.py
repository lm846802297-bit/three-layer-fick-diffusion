from mpmath import iv
from sys import argv
iv.dps=55
I=iv.mpf
pi=iv.pi
T0=I('.0015'); T1=I('.01'); Q=int(argv[1]) if len(argv)>1 else 24; A=12; M=260
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
for j in range(2,6):mom.append((j-1)*mom[j-2]/(2*T0))
U=80/pi*iv.exp(T1*hi**2)*sum(comb(5,j)*hi**(5-j)*mom[j] for j in range(6))
r=iv.exp(-2*pi*A)
alias=2*U*r/(1-r)
# Omitted trapezoid terms, both signs of y.
s=M+1
r=(I(s+1)/s)**5*iv.exp(-T0*(2*s+1))
assert r<1
tail=80/pi*iv.exp(T1*Q**2)*(Q+s)**5*iv.exp(-T0*s*s)/(1-r)
E=(alias+tail).b
print('quadrature error',E,flush=True)
assert E<I('1e-14')
err=I([-E,E])
st=I('.00001'); w=st/2; degree=4
terms=[]; rb=[I(0),I(0),I(0)]
for j in range(M+1):
    q=iv.mpc(Q,j); p=q*q
    z=q/pi*transform(q)*(2 if j else 1)
    amp=iv.exp(p.real*(T1 if p.real>=0 else T0))
    for h in range(3):
        rb[h]+=abs(z)*amp*abs(p)**(h+degree+1)*w**(degree+1)/factorial(degree+1)
    terms.append([z*iv.exp(p*(T0+w)),iv.exp(p*st),[p**z for z in range(degree+3)]])
# Use separate powers of the symmetric cell interval to preserve even-power positivity.
Wp=[I(1)]+[I([0,(w**k).b]) if k%2==0 else I([-(w**k).b,(w**k).b]) for k in range(1,degree+1)]
best=I('-1e20'); low=I('1e20')
for t in range(850):
    vals=[I(0) for h in range(degree+3)]
    for row in terms:
        z,ratio,powers=row
        for m,p in enumerate(powers):vals[m]+=(z*p).real
        row[0]=z*ratio
    jets=[]
    for h in range(3):
        B=rb[h].b
        jets.append(sum(vals[h+k]*Wp[k]/factorial(k) for k in range(degree+1))+I([-B,B])+err)
    f,g,j=jets
    assert f>0,(t,f)
    curv=j/f-(g/f)**2
    assert curv<0,(t,curv)
    best=max(best,curv.b); low=min(low,f.a)
print('a:',a,'b:',b)
print('alias:',alias,'truncation:',tail,'uniform quadrature error:',E)
print('global Taylor errors:',rb)
assert best<I('-2100') and low>I('12.70')
print('850 cells; f above',low,'curvature below',best)
