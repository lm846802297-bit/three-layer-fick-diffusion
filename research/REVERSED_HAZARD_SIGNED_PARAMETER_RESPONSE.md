# Signed transfer response closes the intermediate reversed-hazard window

Keep the three widths (25/66, 4/33, 1/2), diffusivities (25, 1/4, 1), and the release point \(x=1/2\) fixed. Let the left and right Robin rates vary **independently** on \([3.9,4.1]\). The compact-window certificate verifies \(h_X'-h_T'>0\) for \(.0015\le t\le.062\) on the full wall box. Joined to the [analytic origin estimate](REVERSED_HAZARD_FIRST_FACTOR_TAIL_BOX.md) and the [two-pole late-time estimate](REVERSED_HAZARD_COMMON_DENOMINATOR_BRIDGE.md), this gives
\[
h_T'(t)<h_X'(t)\qquad(t>0).
\]
The proof combines analytic bounds and interval arithmetic. It is not a sampled-sign test or a claim for varying layer widths or diffusivities.

## 1. Cancel the common exponential before taking derivatives

Set \(c=\lambda_2,\ a=L_1,\ b=R_1\). If \(N_W\) is the positive normalizing constant of the transferred CDF, use \(U_W=N_W e^{-ct}F_W\); the log curvature is unchanged. On \(q=Q+iy,\ p=q^2\),
\[
 U_W(t)=\frac1\pi\int_{\mathbb R}e^{q^2t}\frac{Y(q)}{D(q)}
           \frac{q^2+\lambda_1}{q^2+r}\,dy ,
\]
where \((Y,r)=(Y_R,b)\) for \(X\) and \((Y_L,a)\) for \(T\). This is the exact inverse Laplace integral: the factor \(q^{-1}\) in the transform cancels against \(dp=2iq\,dy\). Deforming the Bromwich line to the parabola crosses no poles since all singularities in \(p\) are on the negative real axis; transfer decay also suppresses the horizontal connectors. For any real constant \(\beta\), use the centered jets \(V_j=(\partial_t-\beta)^j U_W\). Then
\[
 (\log U_W)''=V_2/V_0-(V_1/V_0)^2.
\]
On each time cell the numerical scale and \(\beta\) are held constant, not differentiated; the scale cancels in both ratios.

## 2. Integrate signed linear responses before applying parameter intervals

Relative to the baseline \((\kappa_L,\kappa_R)=(4,4)\), the amplitude ratio at each contour node has the form
\[
 \frac{(1+y)(1+z)}{(1+d)(1+w)}.
\]
Here \(y\) is the relative change in the selected half trace, \(z=\delta\lambda_1/(p+\lambda_{1,0})\), \(w=\delta r/(p+r_0)\), and \(d=d_1+d_2\), with \(d_1\) linear in the two wall increments and \(d_2\) their product times the normalized mixed coefficient of \(D\). The trace is affine and \(D\) bilinear in the wall rates. The relative amplitude equals \(1+E+R\), exactly, where
\[
 E=y+z-d_1-w,\qquad
 R=\frac{yz-d_2-dw-E(d+w+dw)}{(1+d)(1+w)} .
\]
The script sums each of the four **signed** linear coefficients of \(E\) over the quadrature nodes before multiplying by the interval of its real parameter increment. The two eigenvalue shifts are treated as independent enclosure variables, not as a linearized eigenvalue approximation. Root brackets and corner spectral monotonicity give their enclosures. Only the quadratic residual and other errors are summed in absolute value.

## 3. Certify the whole time cells and the infinite contour

On a time cell \(t=t_m+\xi\) the discrete integrand jet is expanded using \(\exp(p\xi)\) through degree three, with absolute remainder \(\exp(|p||\xi|)(|p||\xi|)^4/4!\). The full interval \(|\xi|\le\Delta t/2\) is evaluated. The residual \(R\) is bounded using \(\exp(|p||\xi|)\) directly.

For the infinite trapezoid rule use \(h=1/2\) and the complex strip \(|\operatorname{Im}y|\le3\). Applying Fourier periodization and shifting the nonzero Fourier modes to the strip edges gives an alias bound: the sum of the two edge \(L^1\) bounds times \(e^{-2\pi\cdot3/h}/(1-e^{-2\pi\cdot3/h})\). The edge integrals are bounded by elementary Gaussian moments. In traveling-wave coordinates the leading coefficient of \(D\) has magnitude at least \(33/16\) after factoring out the leading exponential. The magnitude of every other normalized coefficient is at most one: the interface reflection moduli are smaller than one, and each Robin ratio \((1-\kappa/(\sqrt D q))/(1+\kappa/(\sqrt D q))\) has modulus at most one for \(\Re q>0\). Consequently an absolute lower bound for the coefficient of the scaled denominator on \(\Re q\ge x\) is
\[
 (33/16)\left[2-(1+e^{-2(5/66)x})(1+e^{-2(8/33)x})(1+e^{-x})\right].
\]
Use analogous triangle bounds for the half traces, and \(|q^2+r|\ge x^2\) for \(\Re q\ge x\), \(r>0\). After factoring the Gaussian \(\exp(-ty^2)\), the omitted trapezoid nodes are bounded by polynomial times Gaussian tails with a checked geometric ratio. No high modes are dropped without a tail estimate.

## 4. Scope and reproducibility

Run the companion Python script with mpmath 1.4.1; interval precision is 55 decimal digits. With no arguments the script checks seven closed, contiguous time blocks, each against one full box in the two wall rates. The uniform lower bound on the compact window exceeds 103.85. Lower bounds in the captured stdout are lower bounds for \(h_X'-h_T'\) on the whole indicated block. The origin and late proofs use different estimates and meet the compact proof at its two endpoints.

As a consequence the normalized ratio \(F_T/F_X\) is strictly log-concave, and its log derivative \(h_T-h_X\) decreases strictly to zero. In the same notation, the filtered ratio controlling physical capture has an additional factor proportional to \((h_T-(c-a))/(h_X+(b-c))\). Log-concavity of \(F_T/F_X\) is not by itself a proof that this filtered ratio has a unique peak throughout the wall box, nor a bound on the number of reversals for all three-layer geometries.


## 5. Analytic reduction for the next step

With \(A=c-a>0,\ B=b-c>0,\ \delta=h_T-h_X\), the positive part of the filtered ratio has
\[
r=K\frac{F_T}{F_X}\frac{h_T-A}{h_X+B},
\quad
(\log r)'=\delta+\frac{h_T'}{h_T-A}-\frac{h_X'}{h_X+B}.
\]
Since \(h_T'<h_X'\le0\), \(h_T\) decreases strictly from infinity to zero, so \(r>0\) has one connected component. Also \(\delta\) decreases from infinity to zero (by the unequal optical lengths at the origin), hence it crosses \(A+B=b-a\) exactly once. Writing \(W=h_X+B>0\),
\[
\left(\frac{h_T-A}{W}\right)'=
\frac{\delta'}{W}-\frac{(\delta-A-B)h_X'}{W^2}<0
\quad\text{whenever }\delta\le A+B.
\]
This controls the shape of one factor but does not count critical points of the product. An exact next test, restricted to points of the positive lobe where \((\log r)'=0\), is
\[
(\log r)''=\delta'+\frac{h_T''}{h_T-A}-\frac{h_X''}{h_X+B}
 -\delta\left(\delta-\frac{2h_X'}{h_X+B}\right).
\]
Strict negativity at every such point would imply a unique peak, but is not part of the present certificate.
