# The actual filtered derivative ratio has a unique positive peak

Fixed geometry: widths (25/66,4/33,1/2), diffusivities (25,1/4,1), and kappa=4 at both ends. This is not a statement over a range of parameters.

Let P_L and P_R be the physical cumulative captures, and G_i=P_i+P_i'/lambda_1. The prior ratio reduction gives, for r=G_R'/G_L',

    r(t)=(b/a) M_Z [f_T(t)-(c-a)F_T(t)]/[f_X(t)+(b-c)F_X(t)].

Here a=L_1, b=R_1, c=lambda_2, M_Z>0 is the untilted delay's exponential moment at c, and X,T are the same tilted suffixes as in the reverse-hazard comparison. This is not an exponentially reweighted F_T/F_X. Write N and D_+ for its numerator and denominator without the constant. D_+>0 everywhere. The positive component is (0,t_*), where h_T(t_*)=c-a. Indeed h_T is strictly decreasing by h_T'<h_X'<=0 (F_X is log-concave), with limits +infinity and 0. The peak is unique on this component, with

    0.03438553619384765625 < t_peak < 0.0343856048583984375,
    0.1046212005615234375 < t_* < 0.10462127685546875.

## Exact transform cancellation

With notation from STATIC_TIE_REVERSE_HAZARD_CONTOUR.md, the CDF transforms are positive constant multiples of Y_W(q)(p+lambda_1)/[q D(q)(p+r_W)], with p=s-c, r_T=a and r_X=b. In the transform of N the multiplier is s-(c-a)=p+a. In the transform of D_+ the multiplier is s+b-c=p+b. Thus the half-ground rational factors cancel exactly. The resulting transforms differ only by their constant normalization and the half trace, Y_L vs Y_R:

    U_N(s) = constant_T * Y_L(q) (p+lambda_1)/[q D(q)],
    U_D(s) = constant_X * Y_R(q) (p+lambda_1)/[q D(q)].

No subtraction of nearly equal CDF and density enclosures is needed. The inverse Laplace calculation for the first three jets uses the same pole-free contour, alias estimate and analytic time Taylor remainder as the previous certificate, with polynomial degree six instead of four in the absolute contour bounds. This is because the largest jet weight is (p+lambda_1)(p+c)^2, and the determinant and half-trace bounds are unchanged. The positive normalizations and common exp(ct) factor are not needed for the sign tests.

## Early analytic sign

On 0<t<=0.0015 the earlier saddle-point estimates give h_T-h_X > .037/t^2-412, h_T-(c-a)>.06/t^2 and h_T'>-.128/t^3. For example d_X=7/22, d_T=1/2, x_W=d_W/(2t), and |(log J_W)'|<40/(x_W t), |(log J_W)''|<270/(x_W t^2)+1600/(x_W^2 t^2) for the positive saddle factor J_W. Inserting these bounds in the exact derivative of log r and dropping the nonnegative term -h_X'/(h_X+b-c) yields

    (log r)' > .037/t^2-412-(.128/.06)/t >0.

All rational inequalities for this interval follow by checking its right endpoint after multiplication by the relevant power of t.

## Finite interval sign of the slope

Define Q=N'D_+-ND_+'. On each of 925 closed time cells apply a common gauge beta to the jets, i.e., replace the first jet by f'-beta f and the second by f''-2 beta f'+beta^2 f. A positive common factor was already stripped from all jets. The formulas Q and Q'-2 beta Q are then formed from the gauged jets. The gauge is held fixed on the entire cell and causes no sign change at a zero of Q.

The certificate proves Q>0 through t=.030, Q<0 from .039 through .12, and N>0 throughout [.030,.039]. On every cell in [.030,.039] it proves Q'-2 beta Q<0. Thus at every zero in that window, Q'<0. Endpoints have opposite signs and all zeros cross down strictly, so there is exactly one. The point certificate independently brackets both the peak and the first zero of N. Since h_T is strictly decreasing and that zero occurs before .12, no times after the last finite window belong to the positive component. This proves strict single-peakedness of r on its positive component. No monotonicity claim is made on the negative branch beyond the checked window.

## Consequence for capture crossings at the fixed benchmark

The exact ratio reduction also gives G_L'>0. As t decreases to zero, r tends to zero by the distinct optical exponents. Since (G_R-G_L)'=G_L'(r-1), the level r=1 can be crossed at most twice. In the static tie, G_R-G_L vanishes at t=0 and tends to zero as t tends to infinity. Its derivative is strictly negative both initially and after t_*. It follows that the two crossings of r=1 do exist and that G_R-G_L changes sign exactly once, from negative to positive: it is negative on the initial decreasing segment and positive on the final decreasing segment. The physical cumulative difference obeys (exp(lambda_1 t)(P_R-P_L))'=lambda_1 exp(lambda_1 t)(G_R-G_L). It therefore has at most one positive-time sign reversal. This is an independent upper-bound route for the benchmark, not an upper-bound theorem for other parameters. The actual existence of one reversal also follows from the separately certified late ground-mode sign in FIRST_INVERSION_STATIC_TIE_MIDDLE_CERTIFICATE.md.

Run: python static_tie_physical_ratio.py lo hi Q_X Q_T step mode, using the rows in static_tie_physical_ratio_output.txt. Run the point enclosure with: python static_tie_physical_ratio_points.py .03 .12 16 16 .0001 points. The only Python dependency is mpmath.
