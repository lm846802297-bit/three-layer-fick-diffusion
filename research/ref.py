import mpmath as m
m.mp.dps=35
def capture(t,u,v,z,k,rev,method='dehoog'):
 u,v,z,k,t=map(m.mpf,(u,v,z,k,t))
 def seg(f,p,d,h,q):
  c=m.cosh(d*q);s=m.sinh(d*q)
  return c*f+s*h*p,s*f/h+c*p
 def tr(s):
  q=m.sqrt(s); p=k/q
  pairs=[(z*u,u),((m.mpf('.5')-z)*v,v)]
  if rev:pairs.reverse()
  l,pl=1,p
  for d,h in pairs:l,pl=seg(l,pl,d,h,q)
  r,pr=seg(1,p,m.mpf('.5'),1,q)
  return (k/q)*(r-l)/(s*(l*pr+r*pl))
 return m.invertlaplace(tr,t,method=method)
