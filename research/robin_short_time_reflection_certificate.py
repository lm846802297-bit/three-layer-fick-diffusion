import platform, sys
import mpmath as mp
from mpmath import iv
mp.mp.dps = 90
iv.dps = 70
I = iv.mpf
delta = I(2) / 11
kappa = I(4)
t0 = I('0.0015')
kL1 = I(['3.77267640109840874200470548869422075099999999999999999999999999999999997198','3.77267640109840874200470568869422075100000000000000000000000000000000001776'])
kR1 = I(['4.57785945620680872964282824434745928699999999999999999999999999999999997873','4.57785945620680872964282844434745928700000000000000000000000000000000006073'])
kc = I(['3.97506264914715939438336797080267273199999999999999999999999999999999998388','3.97506264914715939438336817080267273200000000000000000000000000000000002965'])
L1=kL1**2; R1=kR1**2; c=kc**2; b_spectral=R1-L1
L1_upper=I('14.24'); c_upper=I('15.81'); b_upper=I('6.73')
assert L1.b<L1_upper.a and c.b<c_upper.a and b_spectral.b<b_upper.a
a0=4*kappa/5
r_upper=2*t0/delta
rprime_upper=2/delta+4*kappa*t0/delta**2
rsecond_upper=16*kappa/delta**2
lprime_upper=10*t0/delta**2+4*L1_upper*t0**2/delta**2
lsecond_upper=19/delta**2+28*L1_upper*t0/delta**2+4*L1_upper**2*t0**2/delta**2
Ay=2*kappa/delta+4/delta**2
Az=2/delta**2+8*L1_upper*t0/delta**2
EA=Ay+Az
EA2=2*Ay**2+2*Az**2+2*Ay*Az
EB=8*kappa/delta**3+80/delta**4+8*L1_upper/delta**2
jsecond_upper=48*t0/delta**3+48*t0**2*EA/delta**3+8*t0**3*(EA2+EB)/delta**3
jprime_ratio=12*t0/delta**2+4*t0**2*EA/delta**2
jprime_upper=2*t0/delta
j_upper=16*t0**3/delta**3
assert rprime_upper.b<I(12)
assert 2*kappa*t0/delta+2*t0/delta**2<1
assert lprime_upper.b<I(1)
assert t0**2*lsecond_upper.b<(delta**2/2+4*t0).a
assert jprime_ratio.b<I(1)
assert jprime_upper.b<I('0.017')
assert jsecond_upper.b<(12+delta/t0+8/delta).a
H0_lower=1-a0*r_upper-a0*b_upper*j_upper
H0prime_upper=a0*I(12)+b_upper+a0*b_upper*(2*t0/delta)
H0second_t2_upper=t0**2*(a0*rsecond_upper+b_upper*(delta**2/2+4*t0)/t0**2+a0*b_upper*(12+delta/t0+8/delta))
logHsecond_t2_upper=H0second_t2_upper/H0_lower+(H0prime_upper/H0_lower)**2*t0**2
assert H0_lower.a>I('0.94')
assert H0prime_upper.b<I(46)
assert H0second_t2_upper.b<I('0.18')
assert logHsecond_t2_upper.b<I('0.20')
main_coeff=delta**2/2-I('1.7')*t0
assert main_coeff.a>I('0.0139')
main_curvature_upper=-I('0.0139')/t0**3
# The logarithmic first derivative is positive: its t^2-scaled lower
# bound is decreasing in t and remains positive at t0.
logderiv_lower_coeff=delta**2/4-I('1.5')*t0-(c_upper+H0prime_upper/H0_lower)*t0**2
assert logderiv_lower_coeff.a>0
logderiv_coeff=delta**2/4+(c_upper+H0prime_upper/H0_lower)*t0**2
assert logderiv_coeff.b<I('0.009')
# Since (log f0)''<0 and 3/2+(log H0)''*t^2>1.3>0,
# its absolute value is below delta^2/(2 t^3).
assert I('1.5')-I('0.20')>0
logsecond_coeff=delta**2/2
assert logsecond_coeff.b<I('0.019')
a=I(5)/66; b=I(16)/66; d0=delta+2*a
alpha=(d0**2-delta**2)/4
assert alpha.a>I('0.0195')
qmin=d0/(2*t0); sigma0=d0**2/(4*t0**2); contour_a_min=d0**2/(2*t0)
assert contour_a_min.a>I(37)
Mref_bound=(1+b_upper/sigma0)/(1-iv.exp(-qmin))
assert Mref_bound.b<I('1.001')
qfactor=1+c_upper/sigma0
assert qfactor.b<I('1.002')
ct_factor=iv.exp(c_upper*t0)
Jscaled_upper_at_37=iv.sqrt(11*iv.pi)+iv.sqrt(iv.pi)/2*11**(mp.mpf(3)/2)/37+960*11**6/iv.mpf(37)**(mp.mpf(11)/2)
assert Jscaled_upper_at_37.b<I(12)
Cbase=32*d0/(delta*iv.sqrt(2*iv.pi))
eta_coeff=[3*I('1.001')*Cbase*qfactor**n*ct_factor*(d0**2/4)**n for n in range(3)]
assert eta_coeff[0].b<I(73) and eta_coeff[1].b<I('2.1') and eta_coeff[2].b<I('0.057')
assert t0.b<(alpha/7).a
exp_gap=iv.exp(-alpha/t0)
eta0=I(73)*exp_gap
eta1=I('2.1')*t0**-2*exp_gap
eta2=I('0.057')*t0**-4*exp_gap
Hmin=I('0.94')
eps0=eta0/Hmin
L0=I('0.009')/t0**2
Lprime=I('0.019')/t0**3
eps1=eta1/Hmin+L0*eps0
eps2=eta2/Hmin+2*L0*(eta1/Hmin)+(L0**2+Lprime)*eps0
assert eps0.b<I('0.00018') and eps1.b<I(3) and eps2.b<I(50000)
log_reflection=eps2/(1-eps0)+eps1**2/(1-eps0)**2
assert log_reflection.b<I(50000)
full_upper=main_curvature_upper+log_reflection.b
assert main_curvature_upper.b<I('-4100000') and full_upper<I(0)
print('ROBIN_SHORT_TIME_REFLECTION_CERTIFICATE: PASS')
print('delta,d0,alpha',delta,d0,alpha)
print('certified L1,R1,c intervals',L1,R1,c)
print('r bounds',r_upper,rprime_upper,rsecond_upper)
print('l derivative bounds',lprime_upper,lsecond_upper)
print('j bounds',j_upper,jprime_upper,jsecond_upper)
print('H0 lower, |H0prime| upper',H0_lower,H0prime_upper)
print('t0^2 |H0second| upper',H0second_t2_upper)
print('t0^2 |(log H0)second| upper',logHsecond_t2_upper)
print('main log curvature upper',main_curvature_upper)
print('log f0 prime coefficient',logderiv_coeff)
print('log f0 second coefficient',logsecond_coeff)
print('qmin,sigma0,contour a min',qmin,sigma0,contour_a_min)
print('H/N_R factor upper',Mref_bound)
print('qfactor,ct_factor',qfactor,ct_factor)
print('Jscaled upper at a=37',Jscaled_upper_at_37)
print('contour coefficients eta0,eta1,eta2',eta_coeff)
print('exp(-alpha/t0)',exp_gap)
print('epsilon bounds',eps0,eps1,eps2)
print('log reflection curvature upper',log_reflection)
print('full log curvature upper',full_upper)
print('runtime',platform.platform(),'Python',sys.version.split()[0],'mpmath',mp.__version__)
