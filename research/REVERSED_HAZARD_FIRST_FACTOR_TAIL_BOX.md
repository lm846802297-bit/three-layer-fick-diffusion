# A first-delay-factor lower bound and a larger explicit late-time wall box

Fix widths (25/66, 4/33, 1/2), diffusivities (25, 1/4, 1), and release at the midpoint, as in REVERSED_HAZARD_UNIFORM_WALL_BOX.md. In this note the two Robin rates vary independently in [3.9, 4.1]. In the established tilted notation, let \(c=\lambda_2\) and \(\eta=\lambda_3-c\). On the whole box,
\[
 e^{\eta t}\bigl(h_X'(t)-h_T'(t)\bigr)>128,\qquad t\ge0.12.
\]
The separate analytic short-time estimate \(h_X'-h_T'>0.067/t^3\) also holds on this box for \(0<t\le0.0015\). Neither statement covers the intermediate time interval on the larger box. The earlier small-box all-time theorem is unchanged.

## A one-factor certificate for the leading-amplitude gap

Use the normalized CDFs, and write \(1-F_W(t)\sim K_W e^{-\eta t}\). In the ordered suffix setting \(T=X+Z\), the positive delay has the transform
\[
 \mathbb E e^{-s Z}=\prod_{i\ge2}
 \frac{1+s/(L_i-c)}{1+s/(R_i-c)} .
\]
The quadratic eigenvalue growth guarantees convergence at \(s=-\eta\), since \(R_2>\lambda_3\). Comparison of the common-pole residues (or convolution and the exponential-moment bound) gives
\[
 \frac{K_T}{K_X}-1=\mathbb E e^{\eta Z}-1
 \ \ge\frac{\eta(L_2-R_2)}{(L_2-c)(R_2-\lambda_3)}.
\]
Each factor evaluated at \(s=-\eta\) is at least one; keeping only the \(i=2\) factor therefore gives a rigorous lower bound. This avoids subtracting independently enclosed coefficients of the two residues. In particular, any verified lower bound \(\ell\le L_2\), \(\ell>R_2\), may replace \(L_2\) on the right, because \((L-R_2)/(L-c)\) increases in \(L\).

Use the half-phase argument of the small-box note. The left phase bound based on the optical-length gap gives \(R_i<L_i\) for every \(i\ge3\), independently of the positive wall rates here. The interval phase check at \(k=\sqrt{150}\) gives \(L_2>150\), while a root bracket gives \(R_2<110\). The interpolating full roots are enclosed by monotonicity of the Robin quadratic form and absence of root crossings through bracket endpoints. Inserting \(\ell=150\), the interval evaluation gives a relative-gap lower bound greater than 0.2499 throughout the box. The number 0.2499 is a lower bound for \(K_T/K_X-1\), not an estimate of the exact gap.

## Transfer-to-curvature inequality

Suppose for \(t\ge\tau\) that \(F_W(t)=1-K_W e^{-\eta t}+r_W(t)\), where \(|r_W^{(j)}(t)|\le E_{W,j}e^{-\omega t}\), \(j=0,1,2\), and \(\omega>\eta>0\). Set \(y=e^{-\eta\tau}\), \(q=e^{-(\omega-\eta)\tau}\), \(m=1-y(K_X+q E_{X,0})>0\); let \(G\) be any lower bound on \(K_T-K_X\). If \(\eta^2(K_X+G)>q E_{T,2}\), then \(F_T''<0\) throughout the interval, and
\[
 e^{\eta t}(h_X'-h_T')
 \ \ge\ \eta^2G-qE_{T,2}
 -\frac{qE_{X,2}+\eta^2K_X(1-m)}{m}
 -y\left(\frac{\eta K_X+qE_{X,1}}m\right)^2.
\]
Indeed, \(F_T\le1\) and \(F_T''<0\), so \(-F_T''/F_T\ge-F_T''\); drop the nonnegative \(+(F_T'/F_T)^2\) term. Bound the two \(X\)-terms using \(F_X\ge m\). All time-decreasing envelopes can be evaluated at \(\tau\). Put \(G=K_X\) times the one-factor relative lower bound above, retaining the common coefficient rather than subtracting unrelated boxes.

The script uses the established exact transfer contour with central \(\operatorname{Re}p=-110\) and external rays \(q_{\rm contour}=8+iy_{\rm contour}\). The only crossed poles are \(s=0\) and \(s=-\eta\). In fact the split eigenvalues below 110 are \(L_1,R_1,R_2\), and exactly \(L_1,R_1\) lie below the bracket for \(\lambda_3\). Min-max and simplicity therefore permit at most one more full root below 110; the equal determinant signs after the \(\lambda_3\) bracket and at \(\sqrt{110}\) exclude it. The shifted-contour integrated absolute values bound the derivatives of the remainder; \(\omega=110-c\). For \(t\ge\tau=0.12\), the rays are bounded with a Gaussian tail estimate at \(\tau\). The 55-digit interval run gives a strict scaled lower bound exceeding 128.02. No sign or derivative of \(\psi_X\), and no density score of \(Z\), is needed in this late-time proof.

## Origin

The optical lengths do not change on the wall box. In the saddle bound of the earlier theorem, the unmatched Robin factor now deviates by at most \(4.1/x\) and the rational factor by at most \(24/x^2\), where \(x=\operatorname{Re}q_{\rm contour}\ge50\). Relative to the leading exponent, each transfer's reflected part remains bounded by \(7\exp(-2lx)\le1/x\) with \(l=5/66\), so the normalized factor still differs from one by less than \(8/x\): an upper constant is \(100/49+(51/49)(4.1+24/50)<8\). The same Cauchy/saddle derivative estimates and the same optical-gap subtraction in the small-box proof give the short-time bound stated above.

Run `python reverse_hazard_first_factor_tail.py` with mpmath 1.4.1. The default interval is the box in the statement. Outside it, no root identification or conclusion is asserted. Exact direct comparison in the middle interval, rather than proving a single sign change of a different function, is the remaining task for this box.
