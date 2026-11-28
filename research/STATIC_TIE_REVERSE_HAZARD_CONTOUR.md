# Direct reversed-hazard comparison at the fixed static tie

Geometry: widths (25/66,4/33,1/2), diffusivities (25,1/4,1), Robin parameter 4. Let a=L_1, b=R_1, lambda_1 be the full ground eigenvalue, and c=lambda_2, with lambda_1<a<c<b<R_2. Here X and T are the exponentially tilted suffix variables with poles lambda_{i+1}-c for i>=2 and zeros R_i-c and L_i-c, respectively. The reverse hazard is h_W=F_W'/F_W. These are auxiliary distributions, not the physical cumulative capture difference.

Result: h_T'-h_X'<0 for all 0<t<=0.12. On [0.0015,0.12], the computed uniform upper bound is less than -14.51 over 930 closed cells. The sign at sufficiently large time is proved separately without an explicit starting time.

## Exact transforms and the contour

Put p=q^2, l=5/66, n=8/33. Let (Y_L,B_L) and (Y_R,B_R) be the directed Robin solutions at the midpoint, with state (y,D y'/q) taken in the inward direction from either wall. Set D(q)=Y_L B_R+Y_R B_L. The full entire characteristic function in p is q D(q). Its zeros are the negative full eigenvalues. The formulae for the half traces and full determinant are obtained by multiplying the matrices with rows (cosh(dq),r sinh(dq)) and (sinh(dq)/r,cosh(dq)), starting with (1,4/q); here r=1/sqrt(Diffusivity) and d is the optical layer length. Standard separated Sturm-Liouville interlacing through insertion of a midpoint Dirichlet condition identifies the two bracketed full roots. The first split eigenvalues are a and b (the next left eigenvalue lies above b; phase at sqrt(b) lies strictly between pi/2 and 3pi/2). The full roots in (0,a) and (a,b) are lambda_1 and c. The endpoint signs and monotonic phase on the first left-root bracket identify the two first half roots.

The finite optical length gives entire characteristic functions of order 1/2 in p. Genus-zero products over their eigenvalues therefore agree with the transfer-matrix characteristic functions up to constants. Strip the first two full factors and the first half factor and substitute p=s-c. Up to positive constants, the CDF Laplace transforms are

    U_X(s) = Y_R(q)/[q D(q)] * (p+lambda_1)/(p+b),
    U_T(s) = Y_L(q)/[q D(q)] * (p+lambda_1)/(p+a).

The canceled factor p+c is the s in the denominator of the CDF transform. The constants multiply all three jets by the same positive number and drop out of the logarithmic derivative. On q=Q+iy, the three jet integrands, after division by exp(ct), are

    G_j(y,t) = (p+c)^j * [Y(q)/D(q)] * [(p+lambda_1)/(p+r)] * exp(pt)/pi,
    j=0,1,2.

The Bromwich deformation uses fixed-height horizontal connectors in q as in STATIC_TIE_CURVATURE_CONTOUR.md. No full higher eigenvalues need to be evaluated.

The eight exponent terms in D come from the two signs for each of the three layers. The leading coefficient of exp((l+n+1/2)q) is (33/16)(1+4/(5q))(1+4/q). Each other coefficient, relative to this, has absolute value at most one for Re(q)>0: the sign factors are (1+10 sigma_1 sigma_2)/11 and (sigma_3+sigma_2/2)/(3/2), and the Robin ratios have absolute value <=1. Thus for x=Re(q), letting B(x)=(1+exp(-2lx))(1+exp(-2nx))(1+exp(-x)),

    |D(q)| >= (33/16)(2-B(x)) exp((l+n+1/2)x).

The script verifies the constants on each strip x>=Q-A: |Y/D| <= 9 exp(-dx), d_X=l+n=7/22 and d_T=1/2. Each rational ratio is at most 1+21/x^2, because |p+r|>=x^2, and |p+c|<=(1+21/x^2)|q|^2. These bounds give |G_j| <= 10 |q|^4 exp(t(x^2-y^2)-dx) for all three jets. The alias bound uses Poisson summation and the strip shift as in the previous contour certificate. The truncation bound sums the Gaussian-polynomial tail geometrically. All time intervals are covered by point-centered Taylor polynomials with sixth-degree terms and analytic remainder. Before evaluating them, subtract a constant beta times the zero-order jet from the first jet, and replace the second jet by F''-2 beta F'+beta^2 F. This is exact for the curvature formula F''/F-(F'/F)^2, and avoids large interval cancellations. The quadrature and Taylor error estimates are carried through the identical linear combinations.

## Analytic proof all the way to zero

Define A(q) to be e^(dq) Y(q)/D(q) times (p+lambda_1)/(p+r). Normalize by its positive direct-path constant A_infinity. At x>=50, the relative errors in the leading exponent terms of Y and D are each bounded by 1/x. Indeed they are at most 7 exp(-2lx), and 7x exp(-2lx)<1 for x>=50. The remaining Robin factor in A/A_infinity is 1/(1+4/(5q)) for X, or 1/(1+4/q) for T. Its modulus is at most one and its deviation from 1 is at most 4/x. The additional rational factor deviates by at most 21/x^2. Consequently

    |A/A_infinity-1| <= [100/49 + (51/49)(4+21/50)]/x < 8/x.

Use the contour Q=d/(2t), write y=v/sqrt(t) and normalize the Gaussian integral to mass one. Then the exact CDF is a positive constant times

    exp(ct-d^2/(4t)) t^(-1/2) J(t),
    J(t)=E[A(d/(2t)+i v/sqrt(t))/A_infinity],  v has density exp(-v^2)/sqrt(pi).

For 0<t<=0.0015, x=d/(2t)>100 and sqrt(t)/d<1/8 for both choices. The disk of radius x/2 around q stays within x>=50. Cauchy's derivative bounds for A/A_infinity are 32/x^2 and 128/x^3, using the just-proved bound on its deviation from 1. For differentiating under the Gaussian integral,

    |d q/dt| <= (x/t)(1+(sqrt(t)/d)|v|),
    |d^2 q/dt^2| <= (x/t^2)(2+(3 sqrt(t)/(2d))|v|).

The Gaussian moments E|v|<1 and E v^2<1, along with J>=0.92, give |J'/J|<40/(xt) and |J''/J|<270/(xt^2). This provides an explicit bound, not a differentiated formal asymptotic. Since h'=(log F)'', the common 1/(2t^2) term cancels between T and X. The two analytic error bounds, with x=d/(2t), prove

    h_T'-h_X' < -0.067/t^3       for 0<t<=0.0015.

The rational comparison for this constant is included in the script.

## Eventual sign, without a quantitative starting time

Both transforms have their first remaining density pole at -eta, eta=lambda_3-c>0, with strictly positive residue. In fact b<lambda_3<min(R_2,L_2); the strict inequalities follow from Dirichlet splitting and the absence of a shared half root at those boundaries. Equivalently, X begins with a nontrivial atom/exponential mixture of rate eta, and its remaining factors have exponential moments in a neighborhood of eta. The tilted Z suffix has first pole R_2-c>eta. Since T=X+Z independently, the leading density residue of T is that of X multiplied by E exp(eta Z)>1. The transfer estimates at large imaginary s justify shifting the Laplace contour slightly past the first pole, uniformly for the density and its first derivative. Therefore h_T'-h_X' is asymptotic to -eta times the positive difference of the two normalized density residues times exp(-eta t). This gives strict negativity for sufficiently large t but no explicit glue point.

The finite windows and their margins are recorded in the native output. The sign theorem on a bounded window and the eventual sign do not yet exclude a sign change in the remaining gap or prove the original global single-peak statement.

Reproduce a row with: python static_tie_reverse_hazard_contour.py lo hi Q_X Q_T step. The arguments, in this order, are shown in the native output (step is obtained as (hi-lo)/cells). The only external Python dependency is mpmath.
