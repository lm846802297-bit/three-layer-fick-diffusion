# A late-time comparison without a single-crossing hypothesis for psi_X'

Write \(h_X=F_X'/F_X\), \(\psi=f_X'/F_X\), and \(T=X+Z\) independently, with the exact tilted suffix laws. This note sharpens the scalar criterion in research/ROBIN_REVERSE_HAZARD_DERIVATIVE_REDUCTION.md on research/two-parameter-spectral-bridge-2026-11-27. It does not certify the proposed value \(\varepsilon=0.035\), nor does it prove a parameter-global result.

Assume \(h_T\ge h_X>0\), and for \(s\ge\varepsilon>0\) assume
\[
-M\le(\log f_Z)'(s)\le-\rho,\qquad 0<\rho\le M.
\]
For the fixed benchmark these score bounds follow from the certified density log-curvature and its limiting score once \(\varepsilon>z_*\). Let \(t_0\ge2\varepsilon\). A convenient sufficient hypothesis is that \(\psi\) is nondecreasing on \([\varepsilon,\infty)\); nothing about \(\psi'\) is assumed on \((0,\varepsilon)\). More weakly, all that is needed is \(\psi(u)\le\psi(t)\) for \(\varepsilon\le u\le t\) and \(t\ge t_0\). Define \(m=\psi(t_0)\), \(B(v)=m F_X(v)-f_X'(v)\), and \([x]_+=\max(x,0)\). Then
\[
D_{\mathrm{pivot}}=\int_0^\varepsilon
\left(e^{\rho v}[B(v)]_+-e^{Mv}[-B(v)]_+\right)\,dv
\]
satisfies for every \(t\ge t_0\)
\[
h_X'(t)-h_T'(t)\ \ge\
(h_T(t)+h_X(t))(h_T(t)-h_X(t))
+\frac{f_Z(t)}{F_T(t)}D_{\mathrm{pivot}}.
\]

## Proof

The exact conditional expectation identity gives
\[
h_X'-h_T'=(h_T+h_X)(h_T-h_X)+
\frac1{F_T(t)}\int_0^t
[\psi(t)F_X(v)-f_X'(v)]f_Z(t-v)\,dv .
\]
The part \(v\ge\varepsilon\) is nonnegative by the weaker record-value hypothesis. On the other part, replace \(\psi(t)\) by \(\psi(t_0)\), which cannot increase the integrand. For \(v\le\varepsilon\) we have \(t-v\ge\varepsilon\), so the two exponential bounds on the density ratio apply. Multiply the lower weight by the positive part and the upper weight by the negative part; this proves the claim. Both integrals are finite because \(F_X\) and \(f_X'\) are flat at zero. At \(t_0=2\varepsilon\), if the old one-change hypothesis is satisfied, let \(A_0=\psi(2\varepsilon)-\psi(\varepsilon)\ge0\). The pointwise estimate \(m-\psi(v)\ge A_0-\mathbf1_{v<u_0}[\psi(v)-\psi(u_0)]\) shows that \(D_{\mathrm{pivot}}\) is at least the old scalar lower bound. Thus this is a weakening of the old conditions, not merely a reformulation.

For implementation without evaluating the singular quotient near zero, put \(I_\rho=\int_0^\varepsilon e^{\rho v} F_X(v)\,dv\). Two integrations by parts give the exact alternate form
\[
D_{\mathrm{pivot}}
=(m-\rho^2)I_\rho+e^{\rho\varepsilon}
  [\rho F_X(\varepsilon)-f_X(\varepsilon)]
-\int_0^\varepsilon(e^{Mv}-e^{\rho v})
  [f_X'(v)-mF_X(v)]_+\,dv .
\]
The only integral involving a derivative has a continuous integrand at zero. A simple sufficient test for the suffix property is \(f_X'<0\) and \(f_X''\ge0\), since \(\psi'=(f_X''-h_X f_X')/F_X>0\). If \(f_X^{(j)}=(-\eta)^j\eta K e^{-\eta t}+E_j(t)\), \(j=1,2\), with \(|E_j(t)|\le C_j e^{-(\eta+\delta)t}\), those two signs follow for \(t\ge b\) from \(C_j e^{-\delta b}<\eta^{j+1}K\). This reduces a sufficient check to at most a finite interval between \(\varepsilon\) and \(b\) plus a tail estimate. Neither this check at the proposed starting time nor positivity of the pivot integral is asserted here.

## Short-time normalization

For these exact \(X,T\), the CDF Laplace transforms have direct terms proportional to \(s^{-1/2}e^{-d_W\sqrt{s}}\), not \(s^{-1}e^{-d_W\sqrt{s}}\); see STATIC_TIE_REVERSE_HAZARD_CONTOUR.md for the exact \(qD(q)\) denominator. The common short-time CDF power is \(t^{-1/2}\), and thus for each fixed \(t>0\) with \(f_Z\) continuous at \(t\), the conditional CDF near \(u=0\) is \(O_t(u^{3/2} e^{-d_X^2/(4u)})\). These powers do not affect the leading comparison \(h_X'-h_T'=(9/121)t^{-3}+O(t^{-1})\) at the benchmark: the equal \(1/(2t^2)\) terms cancel and the previous differentiated saddle estimates supply the remainder. This is why the assertion is justified by differentiated bounds, not simply by differentiating an undifferentiated big-O expression.
