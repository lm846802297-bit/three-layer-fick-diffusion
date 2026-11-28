# Verified pivot comparison for the fixed kappa=4 benchmark

**Result.** With epsilon=0.035 and t_0=0.070, the pivot integral is greater than 73.9 in the unnormalized transfer convention used below. Moreover f_X'<0 and f_X''>0 throughout [0.035,infinity). Thus the pivot comparison independently proves h_T'<h_X' for all t>=0.070 at this benchmark, using the previously established f_Z score bound. It is an alternative proof of a portion of the already established benchmark theorem, not a new all-time or parameter-global theorem. No hypothesis on the number of sign changes of psi_X' before 0.035 is used.

## Definitions and normalization

The widths, diffusivities, and rates are those of STATIC_TIE_REVERSE_HAZARD_CONTOUR.md. Use p=s-c=q^2. The scripts evaluate the exact transforms
\[
\widetilde{\widehat F}_X(s)=\frac{Y_R(q)}{qD(q)}\frac{p+\lambda_1}{p+b},
\qquad
\widetilde{\widehat f}_Z(s)=\frac{Y_L(q)}{Y_R(q)}\,\frac{p+b}{p+a}.
\]
Tildes suppress positive constants only. Thus psi_X, M and rho are unaffected, and the result for the sign of the integral is unaffected. The reported pivot bound uses this unnormalized F_X. Its omitted normalization is the reciprocal of 0.8452034351... . The bounds for f_X derivatives on the intermediate interval use the same convention. The tail code normalizes F_X to one before testing relative remainder factors.

## Continuous interval certificate

Run:
    python reverse_hazard_pivot_window.py
    python reverse_hazard_pivot_suffix_tail.py

Both runs use 55-digit outward interval arithmetic in mpmath. The first isolates the relevant Robin roots with the same brackets and identification as in the existing transfer proof. It locates R_2 in (10.17^2,10.18^2) by the uniform-right Robin phase. It calculates m=psi_X(0.070), M=-(log f_Z)'(0.035), and rho=R_2-c directly from exact transfers:
- m is between -126.300 and -126.299;
- M is between 89.205 and 89.206;
- rho is between 87.708 and 87.709.

For each positive-time block, invert on q=Q+iy with the strip |Im(y)|<=3. With x=Re(q)>=4, the eight exponential terms give |D| >= (33/16)(2-B(x))exp((d_X+1/2)x), where B(x)=(1+exp(-2lx))(1+exp(-2nx))(1+exp(-x)). The uniform half trace has lower bound (1-exp(-x))exp(x/2)/2, since |(q-4)/(q+4)|<=1. Combining half-trace numerator bounds and |q^2+r|>=x^2 for r>=0, the absolute j-th inversion integrand, including 1/pi, is bounded by 100(1+21/x_min^2)^(j+1)(x+|y|)^(2j+1)exp(-d x + t(x^2-y^2+c)); d is the direct optical distance for the transfer being inverted. The contour is analytic in the strip. Shift the Fourier transform of the integrand to the two strip edges and apply Poisson summation to bound the trapezoid alias by the exponential geometric series. Summing each omitted Gaussian lattice tail by a geometric ratio produces the second explicit error. The floating calculation of the number of quadrature terms is not used to assert any tail tolerance: the actual chosen integer is substituted into the interval error bound. For closed time cells, expand the finite trapezoid sum through degree five at the midpoint; the Taylor remainder for each weight is bounded by exp(|s|*radius)(|s|*radius)^6/6!. Add the continuum quadrature error uniformly over the cell, then use interval rectangle integration of the positive and negative parts. Thus no conclusion is based on time-point sampling.

On [0,0.0015], the short-time saddle estimates of the earlier exact-transfer proof give a uniform bound on |F_X| and |F_X''| using |exp(d_X q)Y_R/D*(p+lambda_1)/(p+b)|<9 at Re(q)=d_X/(2t)>100. With y scaled by sqrt(t), use |q^2+c|<=|q|^2+21 followed by elementary Gaussian moments. Each resulting t^(-k) exp(ct-d_X^2/(4t)), k<=9/2, is increasing on this short interval. This bounds the absolute contribution to the pivot integral by 1.459. The remaining three integral blocks have respective lower bounds -17.300, -22.330, and 115.048. The unrounded interval sum proves D_pivot>73.9.

On the closed cells covering [0.035,0.070], the upper bound on f_X' is -98.01 and the lower bound on f_X'' is 3780.8. On [0.070,0.10] the corresponding bounds are -31.47 and 1181.3. Beyond 0.10 the second script extends the established shifted-contour proof to the third CDF derivative. After removing the pole at -eta, it bounds the remainders relative to the leading terms in the second and third CDF derivatives by 0.060 and 0.854, respectively, uniformly on [0.10,infinity). Because the leading second derivative is negative and the leading third derivative positive, both signs persist. This proves the sufficient suffix monotonicity; the f_Z score bounds follow from the previously certified curvature change and asymptotic score as explained in REVERSED_HAZARD_PIVOT_CRITERION.md.

## Limits of the result

The pivot argument by itself does not address t<0.070. Combining with the previously proved fixed-benchmark short and intermediate pieces gives an alternative assembly of its all-time comparison. The strict numeric margin is not asserted to hold throughout any specified parameter box, and the work does not establish a global bound on physical reversals.
