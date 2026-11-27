import mpmath as m
from scipy.optimize import brentq
m.mp.dps=42
M=m.mpf
def setup(K,d):
 v,u=M(1)/5,M(2)
 R=M('.5')+1/K
 z=(R*(1-d)/(1+d)-1/K-M('.5')*v*v)/(u*u-v*v)
 def tr(k):
  L=(M('.5')-z)*v*k; S=z*u*k; r=k/2
  Y=m.cos(L)+K*v/k*m.sin(L); Q=K*m.cos(L)-k/v*m.sin(L)
  yl=Y*m.cos(S)+u/k*Q*m.sin(S);ql=Q*m.cos(S)-k/u*Y*m.sin(S)
  yr=m.cos(r)+K/k*m.sin(r);qr=K*m.cos(r)-k*m.sin(r)
  return yl*qr+yr*ql,yl,yr
 return tr
def spec(K,d,N):
 tr=setup(K,d)
 def fun(k):return tr(k)[0]
 low=.0001; roots=[]
 while len(roots)<N:
  high=low+.2
  if fun(low)*fun(high)<0:
   r=brentq(lambda x:float(fun(x)),low,high)
   k=m.findroot(fun,(M(r)-M('.00001'),M(r)+M('.00001')))
   _,yl,yr=tr(k)
   c=-2*K*(yr-yl)/(k*m.diff(fun,k))
   roots.append((k*k,c))
  low=high
 return roots
def critical(K):
 (l,c),=spec(K,M(0),1)
 return c
Kstar=m.findroot(critical,(M('7.5'),M('7.6')))
print('kappa star',m.nstr(Kstar,15))
for eps in ['.1','.01','.001','.0001']:
 K=Kstar-M(eps); d=0
 for j in range(4):
  dat=spec(K,d,10); l,c1=dat[0]; c,c2=dat[1]; A=-c2; r=c/l
  x=(c1/(r*A))**(1/(r-1)); t=-m.log(x)/l
  dlead=(r-1)*A*x**r
  def der(t):return sum(l*w*m.exp(-l*t) for l,w in dat)
  tp=m.findroot(der,(t-M('.01'),t+M('.01')))
  dst=sum(w*m.exp(-l*tp) for l,w in dat)
  d=dst
 print('eps',eps,'d_fold',m.nstr(d,9),'time',m.nstr(tp,9),
       'd_fold/d_leading',m.nstr(dst/dlead,12), 'time correction',m.nstr(tp-t,6))
