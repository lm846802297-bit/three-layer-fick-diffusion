"""High-precision illustration (not used for any proof)."""
from late_triple import eq, make
import numpy as np
from mpmath import mp as M
M.dps=80
f=M.mpf
h0,K0,l1,l2=M.findroot(lambda *a:eq(a,M),('0.204488687014','7.617134583573','7.184051634522','26.037661337960'))
def modes(q,h,K,d):
    z=f('.5')-h
    rho=(f('.5')+1/K)*(1-d)/(1+d)
    u=M.sqrt((rho-1/K-h/4)/z)
    y,p,n=1,K/q,0
    for j,(w,r) in enumerate(((h,f('.5')),(z,u),(f('.5'),f(1)))):
        a,b,x=y,r*p,w*r*q
        C,S=M.cos(x),M.sin(x)
        n+=w*(a*a+b*b)/2+((a*a-b*b)*S*C/2+a*b*S*S)/(r*q)
        y,p=a*C+b*S,(b*C-a*S)/r
        if j==1:mid=y
    return p+K*y/q,K*mid*(1-y)/(n*q*q),mid*(1+y)/n
def eig(q,h,K,d):
    return M.findroot(lambda k:modes(k,h,K,d)[0],(q-f('.01'),q+f('.01')))
k3=eig(f('7.4925'),h0,K0,0)
l3=k3*k3
c3=modes(k3,h0,K0,0)[1]
offset=list(map(f,('.02','.04','.06')))
alpha=M.lu_solve(M.matrix([[1,-M.exp(-l1*u),-M.exp(-l2*u)] for u in offset]),
                M.matrix([c3*M.exp(-l3*u) for u in offset]))
w=list(map(f,('0','.03','.05','.08')))
g=lambda t:alpha[0]-alpha[1]*M.exp(-l1*t)-alpha[2]*M.exp(-l2*t)-c3*M.exp(-l3*t)
p=lambda x:M.nstr(x,12)
print('lambda:',*(p(x) for x in (l1,l2,l3)))
print('alpha:',*(p(x) for x in alpha))
print('g witnesses:',*(p(g(t)) for t in w))
S1=modes(M.sqrt(l1),h0,K0,0)[2]*K0/l1
print('survival principal coefficient',p(S1))
for T in map(f,('.25','.5','1','2')):
    d=alpha[0]*M.exp(-l3*T)
    ca,cb=[alpha[i]*M.exp(-(l3-l)*T) for i,l in [(1,l1),(2,l2)]]
    def res(h,K,q1,q2):
        a=modes(q1,h,K,d); b=modes(q2,h,K,d)
        return a[0],b[0],a[1]-ca,b[1]-cb
    h,K,q1,q2=M.findroot(res,(h0,K0,M.sqrt(l1),M.sqrt(l2)),tol=f('1e-140'))
    z=f('.5')-h
    r=M.sqrt(((f('.5')+1/K)*(1-d)/(1+d)-1/K-h/4)/z)
    gs=np.sqrt(make(float(K),z=float(z),u=float(r),v=.5,n=24)[1])
    q=[q1,q2]+[eig(f(str(x)),h,K,d) for x in gs[2:]]
    assert len(q)==24 and all(b>a for a,b in zip(q,q[1:]))
    cs=[modes(k,h,K,d)[1] for k in q]
    def delta(t):
        return d-sum(c*M.exp(-k*k*t) for k,c in zip(q,cs))
    rr=[M.findroot(delta,(T+v-f('.002'),T+v+f('.002')),tol=f('1e-130')) for v in offset]
    print('T',p(T),'offsets',*[p(r-T) for r in rr],
          'Delta(T+.03)',p(delta(T+f('.03'))),
          '|Delta|/Survival (leading S)',p(abs(delta(T+f('.03')))/(S1*M.exp(-l1*(T+f('.03'))))))
