# Analytical stability of the derivative peak and the capture crossing

## Statement and scope

Start from widths (25/66,4/33,1/2), diffusivities (25,1/4,1), release at 1/2, and Robin reaction coefficients (4,4). Allow nearby positive layer diffusivities, a small displacement of the internal interface in the left half (so its two widths still sum to 1/2), and independent perturbations of the two Robin coefficients. The release and the right-hand interface remain at 1/2. There is an open neighborhood of the baseline with the following properties.

1. The physical filtered derivative ratio r=G_R'/G_L' has a single positive component and a unique strict maximum there. Both the maximum and the end of the positive component persist as simple critical/zero points.
2. For the convention Delta=P_R-P_L, let S=Delta(infinity). There is exactly one simple capture reversal if S>=0, and exactly two simple reversals if S<0. The second reversal in the latter case escapes to infinity as S tends to zero.
3. If S=-epsilon<0, the escaping reversal satisfies the two-term asymptotic below, with coefficients evaluated at the perturbed parameters.

The propagation from the baseline is analytic. It uses the already established strict fixed-case signs and simple zero/critical points as input, some of which were computer-assisted. This is not a fresh purely analytic proof of that input, nor an explicit neighborhood or a global family theorem. The one-sided late-root mechanism for common Robin coefficients already appears in FIRST_INVERSION_POSITIVE_STATIC_FOLD_LAW.md (with the opposite convention for Delta). The present argument gives the open-neighborhood stability of the actual derivative peak, allows independent wall reaction coefficients, and gives a higher-order escaping-root expansion.

## Two uniform endpoint controls

For a directed half solution initialized by (u,D u'/q)=(1,kappa_i/q) at wall i, let Y_i and B_i be its two components at the release point, and write Dcal=Y_L B_R+Y_R B_L. The (unfiltered) right and left flux transforms are kappa_R Y_L/(q Dcal) and kappa_L Y_R/(q Dcal), respectively, with q=sqrt(s). These follow from Green-function matching and the unit flux jump at the release point. The two optical lengths are d_L=sum_{left} widths/sqrt(D) and d_R=sum_{right} widths/sqrt(D). At the baseline d_L=7/22<d_R=1/2, and this ordering is strict nearby.

Fix a sufficiently small closed parameter box temporarily, only to obtain uniform bounds. Each layer transfer is a sum of two exponentials times constant matrices; at the baseline and nearby the leading coefficient of Dcal (all propagation signs positive) is strictly positive and bounded below as Re(q) grows. In each normalized transfer ratio H_i(q)=exp(d_i q)Y_other(q)/Dcal(q) the leading limit H_i(infinity)>0, and after factoring the leading exponential the other terms have bounds O(exp(-delta Re(q))) for a common delta>0. The boundary factors differ from their limits by O(1/|q|). Consequently H_i/H_i(infinity)=1+O(1/x), uniformly for x=Re(q) sufficiently large and all Im(q). Cauchy's estimate on the disk of radius x/2 gives a uniform O(x^-2) bound on its q derivative.

Invert on the contour q=d_i/(2t)+iv/sqrt(t). This is the same deformation as in STATIC_TIE_REVERSE_HAZARD_CONTOUR.md, valid uniformly since the optical gap stays positive and all remaining terms have bounded exponential coefficients. The Jacobian cancels the q in the flux transform. Thus the flux f_i=P_i' equals a positive smooth parameter-dependent constant times

    t^(-1/2) exp(-d_i^2/(4t)) J_i(t),    J_i(t)=1+O(t).

The last bound follows by integrating against the fixed Gaussian exp(-v^2), since x=d_i/(2t). In G_i'=f_i+f_i'/lambda_1 the transform is multiplied by 1+q^2/lambda_1, with no initial term since f_i(0+)=0. Factor x^2/lambda_1 on the same contour. It follows that G_i' is a positive constant times

    t^(-5/2) exp(-d_i^2/(4t)) K_i(t).

The Gaussian integrand for the normalized K_i is H_i/H_i(infinity) times [(1+2iv sqrt(t)/d_i)^2+4 lambda_1 t^2/d_i^2]. Therefore K_i=1+O(sqrt(t)), and K_i'=O(t^-1/2). The derivative bound follows by differentiation under the Gaussian: H_i' is O(x^-2), q_t is bounded by x/t+|v|/(2t^(3/2)), and differentiating the displayed polynomial costs at most O(t^-1/2)(1+v^2). All constants are uniform on the small box. In particular G_R' and G_L' are both positive for uniformly small t, r tends to 0, and

    (log r)' = (d_R^2-d_L^2)/(4t^2) + O(t^-1/2)>0.

The unfiltered flux ratio f_R/f_L also tends uniformly to 0, so Delta'<0 and Delta<0 initially.

For late time use the normalized eigenfunction expansion. Min-max comparison with the Neumann Laplacian bounds lambda_n below by a uniform positive constant times (n-1)^2. The elementary one-dimensional Sobolev bound on eigenfunctions (from normalization and the quadratic form) bounds their wall and release values by a fixed polynomial in lambda_n. Thus all time derivatives needed on a positive compact time interval are given by locally uniformly convergent series. For each fixed n, the eigenvalues and spectral coefficients depend smoothly on the parameters: use the transfer matrix in each layer and the simplicity of the separated eigenvalues, or pull the pieces back to fixed intervals. These facts also give uniform late estimates with the next-mode exponential rate, by factoring it from the absolutely summable tail at one fixed positive time. In the individual G_i' the ground mode is canceled identically. The second-mode amplitudes are nonzero and have opposite signs at the baseline: positive on the left and negative on the right. This follows from the tilted CDF formulas in STATIC_TIE_PHYSICAL_RATIO_PEAK.md and their positive normalizations, since b-c>0 and a-c<0. The sign of each amplitude and a positive gap to mode 3 persist. Hence G_L'>0 and G_R'<0 for all uniformly late t.

For Delta itself the ground mode does not cancel. Write

    Delta = S + A exp(-lambda_1 t) + B exp(-lambda_2 t) + O(exp(-lambda_3 t)).

The first coefficient A is strictly positive at the baseline in the convention P_R-P_L. This is the ground-mode contrast in FIRST_INVERSION_STATIC_TIE_MIDDLE_CERTIFICATE.md after reversing its Delta sign; the positive ground state has larger left-wall flux. A stays bounded below and the expansion including its first time derivative is uniform. Thus Delta is strictly decreasing for all uniformly late t, on both sides of the tie.

## Finite-time continuation and exact count

Use a common short cutoff and late cutoff from the two endpoint arguments. At the baseline the denominator G_L' is positive for all times. The ratio strictly increases then strictly decreases on its positive component; its maximum is nondegenerate, the zero ending its positive component is simple, and the ratio is strictly negative afterwards. These are the signs of the fixed-case certificate and the global reverse-hazard result. On the remaining compact window, continuity in the parameters preserves all strict signs away from small disjoint intervals around those two points. In those intervals the corresponding derivative (r'' at the peak, r' at the zero) keeps its strict sign by compact C^2 continuity. Endpoint signs and the implicit function theorem then give exactly one maximum and exactly one ending zero throughout a sufficiently small parameter neighborhood. The endpoint controls exclude any new positive component at zero or infinity.

For Delta, the static-tie baseline has one negative-to-positive simple reversal. Simplicity also follows from its filtered representation: at its zero the weighted derivative is strictly positive, because G_R-G_L has already changed from negative to positive. The same compact continuation gives exactly one finite-time zero in the common window, with Delta positive at its late endpoint. The uniform decreasing late tail then has no finite crossing if S>=0 (including S=0), and exactly one simple downward crossing if S<0. Negative early sign is uniform. This argument does not presume that the physical capture difference is itself a ratio.

The sign of S has a closed static expression. Let R_i=1/kappa_i+integral_half_i (1/D) dx be the diffusive resistance plus wall resistance, so

    P_R(infinity)=R_L/(R_L+R_R),    S=(R_L-R_R)/(R_L+R_R).

Indeed, the steady splitting function has constant D u', and its Robin boundary conditions are D u'(0)=kappa_L u(0) and D u'(1)=kappa_R (1-u(1)); solving these gives the formula at the release point. In the baseline geometry, perturb only the right wall coefficient. Then S=(1/4-1/kappa_R)/(5/4+1/kappa_R), so the second reversal exists for kappa_R<4 near 4; it disappears to infinity at the tie, not by collision with the first reversal.

## Two-term escape law

Put epsilon=-S>0, x=exp(-lambda_1 t), p=lambda_2/lambda_1>1, q=lambda_3/lambda_1>p, and x_0=epsilon/A. At the late root the exact equation has the form

    -epsilon+A x+B x^p+O(x^q)=0.

The smooth, uniform spectral tail implies x/x_0 -> 1. Taylor expanding x^p at x_0 with p in a compact positive interval gives

    x=x_0-(B/A)x_0^p+O(x_0^min(2p-1,q)),
    t_late=(1/lambda_1) log(A/epsilon)
             +(B/(lambda_1 A)) x_0^(p-1)
             +O(x_0^min(2p-2,q-1)).

All coefficients and exponents are those of the perturbed system; local uniformity follows from positive lower bounds on A and both spectral gaps. This avoids an invalid substitution of the baseline exponent along paths approaching the tie surface tangentially. For the transverse one-wall path, epsilon=(4-kappa_R)/24+O((4-kappa_R)^2), hence the leading law can equivalently use baseline data:

    t_late=lambda_1,0^(-1) log[24 A_0/(4-kappa_R)]+o(1).

The neighborhood size is existential. No additional parameter grid or interval test is used here.
