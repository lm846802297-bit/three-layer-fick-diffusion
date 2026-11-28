# All-time reversed-hazard derivative comparison: fixed static-tie benchmark

For the fixed geometry in STATIC_TIE_REVERSE_HAZARD_CONTOUR.md, the reversed hazards h_W=(log F_W)' satisfy

    h_T'(t) < h_X'(t)   for every t>0.

The antecedent certificate covers (0,0.12]. The following proof covers [0.12,infinity), including all the modes beyond the common slowest mode. This is a fixed-parameter, computer-assisted theorem, not a uniform statement in the geometry.

## The shifted contour

Use the exact CDF Laplace transforms from the finite-time proof: U_W(s)=Y_W(q)/(q D(q)) (p+lambda_1)/(p+r_W), p=s-c=q^2, with r_X=b and r_T=a. The choice Y_X=Y_R and Y_T=Y_L is unchanged. Let K_W=Res_{s=0} U_W(s)>0; then U_W/K_W is the normalized CDF transform. The two full roots below the split threshold b are removed from the density transforms, and the half ground roots at a,b cancel their respective denominators. The first remaining full pole is lambda_3, with eta=lambda_3-c>0. We move left of that pole to Re(p)=-B, B=110, without crossing lambda_4.

Here is an explicit check on the pole count. In k>0 coordinates there is a full root below sqrt(a), one in (sqrt(a),sqrt(b)), and one between 7.29 and 7.30. Inserting a Dirichlet condition at the midpoint gives interlacing with the sorted two half spectra. There are exactly two split eigenvalues below 7.30^2, and exactly three below B, namely a,b,R_2; this is certified by the increasing left phase and the exact right root equation. In particular, the left phase at sqrt(B) is between pi/2 and 3pi/2, and sqrt(B) is between the second and third right roots. Thus at most one further full root can be below B. The signed determinant immediately to the right of the third-root bracket and at sqrt(B) has the same strictly positive sign; since full roots are simple, there is no such extra root. All comparisons required for these brackets are interval-checked. The sign convention for the determinant is Im(D(ik)).

For the upper half of the new contour use q=x+i sqrt(B+x^2) for 0<=x<=8, then q=8+iy for y>=y_0=sqrt(B+64). The other half is the conjugate. The straight central contour is Re(p)=-B. On the bent tails Re(p)=-B-(y^2-y_0^2). Moving the vertical p-line first and then bending these tails is valid: at large imaginary p the transfer ratio decays as exp(-d Re(sqrt(p))) times a polynomial, and the horizontal q connectors on bending have Gaussian decay for Re(q)<=Im(q)/2 and exponential optical decay elsewhere. For j=0,1,2 the integrand after the Jacobian contains (p+c)^j Y_W/D (p+lambda_1)/(p+r_W), divided by pi K_W. This is the same analytic integrand as before, with the explicit residues added.

## Certified amplitude and remainder

Let A_W=-Res_{s=-eta} U_W/K_W. Numerically enclosed values are

    eta = 37.42744995205...,
    A_X = 1.12090953989...,
    A_T = 2.01864504698....

These are survival amplitudes, not density amplitudes. The spectral derivative in the residue is exact: if H(k)=Im D(ik), then (q D(q))'_p at p=-k^2 is H'(k)/2 at a full root. The script differentiates the real transfer formula analytically.

For j=0,1,2 write

    F_W^(j)(t) = 1_{j=0} - A_W (-eta)^j exp(-eta*t) + E_Wj(t),
    |E_Wj(t)| <= M_Wj exp(-(B-c)*t),   t>=t0=0.12.

The central part of each M_Wj uses a modulus bound on every one of 2048 closed contour intervals (integrated using its width), not point quadrature. For the infinite vertical tails use the explicit eight-term determinant lower bound at Re(q)=8 and the separate leading-exponent bounds on the two half traces. After y=y_0+u, retain exp(-2*y_0*t0*u) and discard the additional exp(-t0*u^2). Integrating the remaining positive polynomial is exact via its elementary exponential moments. This bounds the entire omitted part with no high mode enumeration.

Finally write s0=exp(-eta*t0), e=exp(-(B-lambda_3)*t0). Positivity and F_W<=1 give

    F_W in [1-(A_W+M_W0*e)*s0, 1],
    exp(eta*t) f_W' in -eta^2 A_W + [-1,1]*M_W2*e,
    exp(eta*t) f_X <= eta A_X+M_X1*e.

With nonnegative densities the square of the last upper bound bounds the square of the scaled density. Substitute these enclosures into h_T'-h_X' = f_T'/F_T-f_X'/F_X-(f_T/F_T)^2+(f_X/F_X)^2, drop the nonpositive T square, and scale by exp(eta*t). The script verifies

    exp(eta*t)*(h_T'-h_X') < -1185    for all t>=0.12.

This overlaps the finite-time proof at its exact endpoint. Consequently F_T/F_X is strictly log-concave on (0,infinity). Its log derivative decreases from +infinity (the difference of the two optical short-time exponents) to 0 (the common spectral tail); hence the unweighted ratio increases, and every ratio exp(-d*t)*F_T/F_X with fixed d>0 has exactly one interior maximum. This does not identify the original physics observable with such a ratio; that identification or any extra factors must be checked separately.

Reproduce with: python static_tie_reverse_hazard_tail.py (dependency: mpmath).
