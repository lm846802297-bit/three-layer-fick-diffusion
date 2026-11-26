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
if __name__ == '__main__':
    m.mp.dps=45
    for t in ('0.01','0.05','0.13','2'):
        for fast_at_wall in (True,False):
            vals=[capture(t,'1.1','.04','.4132','6.5',fast_at_wall,alg)
                  for alg in ('dehoog','cohen')]
            assert abs(vals[0]-vals[1]) < m.mpf('1e-29')
            print(t,'fast at wall' if fast_at_wall else 'slow at wall',
                  *(m.nstr(x,17) for x in vals))
