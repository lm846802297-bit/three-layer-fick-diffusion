# Uniform reverse-hazard comparison near the two-wall benchmark

**Parameters.** Keep the three widths (25/66,4/33,1/2), the diffusivities (25,1/4,1), and the release point at the middle interface fixed. Let the Robin reaction rates at the left and right walls vary independently:
\[
(\kappa_L,\kappa_R)\in [3.99999,4.00001]^2.
\]
**The result is \(h_T'(t)<h_X'(t)\) for every \(t>0\) throughout the box.** The comparison concerns the two tilted suffix CDFs, with \(c=\lambda_2\) for each geometry, and \(h_W=(\log F_W)'\).

## Why the positive-delay decomposition persists on the whole box

For \(i\ge2\), denote the \(i\)-th left and right half eigenvalues by \(L_i,R_i\). Their roots are those of the half problems with Dirichlet condition at the release point. The right \(i\)-th wave number lies between \((2i-1)\pi\) and \(2i\pi\), by the Robin phase equation. The unwrapped left phase is
\[
\theta_L(k)=k(8/33)+
\operatorname{arg}_{\rm lift}\!\bigl(\cos v+10\,\mathrm i\sin v\bigr),
\quad v=k(5/66)-\arctan(\kappa_L/(5k)).
\]
The lift is chosen to agree with zero at \(v=0\). Its difference from \(v\) is strictly less than \(\pi/2\) in absolute value, since the scalar product of \((\cos v,\sin v)\) and \((\cos v,10\sin v)\) is strictly positive. Consequently \(\theta_L(\sqrt{R_i})<(7/11)i\pi+\pi/2 < (i-1/2)\pi\) for every \(i\ge3\). The left half phase is strictly increasing in \(k\); its \(i\)-th crossing is of \((i-1/2)\pi\). For \(i=2\), the interval phase check at \(\sqrt{110}\) gives \(\theta_L(\sqrt{110})<3\pi/2\) and \(\sqrt{R_2}<\sqrt{110}\). Hence \(R_i<L_i\) for every \(i\ge2\). The first-mode inequalities \(\lambda_1<a<c<b<\lambda_3\), \(a=L_1\), \(b=R_1\), are retained.

For \(i\ge2\), there are at least \(i+1\) split eigenvalues at or below \(R_i\), counting \(L_1\) and the first \(i\) right modes. The min-max principle therefore gives \(\lambda_{i+1}\le R_i\). The tilted product for \(F_X\) is the law of a sum of independent atom/exponential mixtures. Each such law has log-concave CDF: on \(t\ge0\) the CDF is \(1-qe^{-\omega t}\), \(0\le q\le1\), and is zero on \(t<0\). This property survives summation: write each variable as the generalized quantile of an independent uniform, equivalently as a nondecreasing convex function of a negative exponential, by inversion of its concave log CDF. The sum sublevel set is convex; marginal log-concavity proves the assertion for finite sums, then pass to the limit (finite sum of expectations by quadratic eigenvalue growth). The ordered delay factors also give an independent nonnegative \(Z\) with \(T=X+Z\). The decreasing reverse hazard of \(X\) implies \(h_T\ge h_X\) by conditioning, as required for the pivot criterion.

## Verification on a continuum of wall parameters

The interval code substitutes the box interval at every occurrence of a wall reaction rate. Even where the same interval object is reused for the two walls, this encloses every independent pair, not just the diagonal. For roots, the endpoint brackets do not lose their sign changes anywhere in the box. Each indexed eigenvalue increases in each reaction rate by its quadratic-form min-max formula, so the two equal-rate corners bound all intermediate roots. Fixed-saddle time blocks use the same strip and Gaussian alias/tail method as the fixed-benchmark certificate. For compact curvature comparisons, the CDF jets are centered before enclosing the time variation: for a fixed cell-dependent \(\beta\), use \(F,\,F'-\beta F,\,F''-2\beta F'+\beta^2F\). Their curvature formula equals the uncentered one.

For all parameters in the box:
- On [.035,.070] the score combinations \(f_Z'+85f_Z\) and \(f_Z'+95f_Z\) are respectively negative and positive by compact contour bounds on the combinations themselves.
- On [.070,infinity) the same signs follow by shifting the delay-density contour past its first remaining pole \(-\rho\), with \(\rho=R_2-c\). In \(p=s-c\) the new contour crosses only the second right root (the first cancels). Use \(\operatorname{Re}p=-220\) in the central segment and \(q=8+iy\) on the external rays. The scaled error / leading margin ratios are <.053 and <.019 for the two score combinations and decrease thereafter. These joint inequalities also imply \(f_Z>0\) on [.035,infinity). No parameter-uniform uniqueness of the curvature zero of the delay density is assumed.
- The density derivatives \(f_X'<0<f_X''\) hold on [.035,.10] by full time cells and on [.10,infinity) by the first-pole bound, with relative errors <.060 and <.853.
- Thus the fixed weights 85 and 95 may be used as \(\rho\) and \(M\) in the interval form of the pivot criterion (85 is a score lower rate bound, not the actual pole). For \(\varepsilon=.035,t_0=.070\), the unnormalized pivot lower bound exceeds 64.14. This proves \(h_T'<h_X'\) for \(t\ge.070\).

## Analytic control all the way to zero

The previous short-time proof is uniform here. In the saddle representation
\[
F_W(t)=C_W e^{ct-d_W^2/(4t)}t^{-1/2}J_W(t),
\quad d_X=7/22,\ d_T=1/2,
\]
the normalized analytic transfer factor obeys \(|A(q)/A_\infty-1|<8/\operatorname{Re}q\) for \(\operatorname{Re}q\ge50\). Indeed the nonleading exponent terms relative to the leading ones have total absolute value at most \(7e^{-2l\operatorname{Re}q}\) (\(l=5/66\)); the two Robin reflection moduli are at most one as the wall rates are positive, and the unmatched Robin factor differs from one by at most \(4.01/\operatorname{Re}q\). The remaining rational factor differs by at most \(21/(\operatorname{Re}q)^2\). Thus the same Gaussian and Cauchy derivative estimates apply simultaneously for every pair of wall rates. In particular, with \(x=d_W/(2t)>100\) for \(t\le0.0015\),
\[
|J_W'/J_W|<40/(xt),\qquad
|J_W''/J_W|<270/(xt^2).
\]
The common \(1/(2t^2)\) terms cancel in the difference of log curvatures. Subtracting the uniform two error bounds from \((d_T^2-d_X^2)/(2t^3)\) gives \(h_X'-h_T'>0.067/t^3\) for \(0<t\le0.0015\).

## Join and conclusion

On [0.0015,0.070] the largest outward upper bound over all closed time cells is less than -137.52 (for h_T'-h_X'). Together with the analytic origin estimate and the uniform pivot estimate, this covers all positive times without gaps. It strengthens the earlier existence of an open neighborhood by giving a concrete product box for independent wall parameters. It is a proof using analytic estimates and outward interval arithmetic, not a grid of parameter sample points.

For reproduction with Python 3.14.7 and mpmath 1.4.1, run the two scripts without command-line arguments. The default box is the one stated above. A candidate run for a larger box must still satisfy the analytic hypotheses, all interval assertions, and the spectral root identifications; no larger box is asserted here. Neither changes to layer thicknesses/diffusivities nor a general bound on numbers of capture reversals is asserted.
