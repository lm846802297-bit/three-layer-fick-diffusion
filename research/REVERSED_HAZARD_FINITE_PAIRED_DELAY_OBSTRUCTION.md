# The finite paired Robin product reverses the short-time curvature comparison

## Setting and result

Write \(h_W=(\log F_W)'\) for the reverse hazard, and use the tilted suffix notation \(T=X+Z\) with independent summands. Keep the exact \(X\), and truncate **only the paired delay factors** of \(Z\) at index \(N\). In the first-inversion regime with ordered suffix \(L_i\ge R_i>c\), each tilted factor is
\[
\frac{1+s/(L_i-c)}{1+s/(R_i-c)}
=p_i+(1-p_i)\frac{\rho_i}{s+\rho_i},\qquad
\rho_i=R_i-c,\quad p_i=\frac{R_i-c}{L_i-c}.
\]
This follows by evaluating the untilted factor at \(s-c\) and normalizing its value at \(s=0\). Let \(Z_N\) be the sum of these independent mixture delays for \(2\le i\le N\), and \(T_N=X+Z_N\).

Put \(A=d_X^2/4>0\). For each **fixed** \(N\) with at least one strict pair,
\[
\boxed{\displaystyle h_{T_N}'(t)-h_X'(t)
=\frac{2}{A}\sum_{i=2}^{N}(L_i-R_i)+O_N(t)\qquad(t\downarrow0).}
\]
Thus every nontrivial finite paired product gives strictly the **opposite** short-time sign from the full product, for which \(d_T>d_X\) yields
\[
h_T'(t)-h_X'(t)=-\frac{d_T^2-d_X^2}{2t^3}+O(t^{-1})
\]
by the differentiated short-time transfer estimate. This is a limitation of this truncation procedure, not a counterexample for the diffusion and not a statement about all spectral truncation procedures.

## Proof

The law of \(Z_N\) has an atom \(p_N^{(0)}=\prod_{i=2}^N p_i>0\) at zero and an absolutely continuous part with density \(g_N\), smooth from the right at zero. Convolutions involving two or more exponential densities vanish at zero. Consequently,
\[
\frac{g_N(0+)}{p_N^{(0)}}=\sum_{i=2}^N\frac{1-p_i}{p_i}\rho_i
=\sum_{i=2}^{N}(L_i-R_i).
\]
Assume the previously established short-time form \(F_X(t)=C e^{ct-A/t}t^{-1/2}J(t)\), with \(J=1+O(t),\ J'=O(1),\ J''=O(t^{-1})\); see REVERSED_HAZARD_PSI_OBSTRUCTION_AND_LOCAL_STABILITY.md. In the convolution ratio, use \(z=t^2y/(1+ty)\), so \(u=t-z=t/(1+ty)\) and \(1/u-1/t=y\). The formula is exact:
\[
\frac{F_{T_N}(t)}{F_X(t)}
=p_N^{(0)}+t^2 M_N(t),\quad
M_N(t)=\int_0^\infty
e^{-Ay}(1+ty)^{-3/2}e^{-cz}
\frac{J(u)}{J(t)}g_N(z)\,dy.
\]
For small \(t>0\), the first \(t\)-derivative of the integrand is dominated by a fixed polynomial in \(y\) times \(e^{-Ay}\); for the second derivative the same holds after multiplication by \(t\). In particular \(u'=(1+ty)^{-2}\), so the potentially singular term obeys \(|J''(u)|(u')^2\le C/[t(1+ty)^3]\); all other derivatives are bounded by polynomials in \(y\) since \(u,z\in[0,t]\), \(J(t)\) stays away from zero, and \(g_N\) has bounded derivatives there. Hence \(M_N=g_N(0+)/A+O_N(t)\), \(M_N'=O_N(1)\), and \(M_N''=O_N(1/t)\). Differentiating the logarithm of \(p_N^{(0)}+t^2M_N(t)\) twice now proves the asserted formula. In particular, this calculation does not differentiate an uncontrolled remainder.

## Consequences for a proof strategy

1. The mechanism is the artificial zero-delay channel. For fixed \(N\), the atom has the same optical exponent as \(X\); the absolutely continuous part first contributes at relative order \(t^2\). The full delay law has no such shortcut when \(d_T>d_X\). For piecewise constant positive media, Dirichlet/Neumann bracketing on the finitely many layers gives \(L_i\sim (\pi i/d_X)^2,\ R_i\sim (\pi i/d_T)^2\), so \(p_N^{(0)}\) decays exponentially and the coefficient \(\sum_{i=2}^N(L_i-R_i)\) grows cubically. These estimates do not identify the location of the spurious transition.
2. There are actual spurious curvature zeros for large truncations. The expected sum of the delays is finite since \(\sum_i1/\rho_i<\infty\), so \(Z_N\to Z\) almost surely under the independent-sum coupling. For each \(t>0\), the convolution and its first two derivatives converge by dominated convergence: the zero extension of each \(F_X^{(j)}\), \(j=0,1,2\), is continuous and bounded on \((-\infty,t]\) by the short-time estimate. Therefore the curvature difference converges pointwise for \(t>0\). At any sufficiently small fixed positive \(t\) the full sign is negative, hence eventually so is the finite sign there, whereas for each fixed nontrivial truncation the sign is positive sufficiently near zero. Continuity implies a sign-changing transition between the two. The first positive zero tends to zero: the same argument works with an arbitrarily small fixed right endpoint. No rate is established.
3. For the unsolved parameter region, this rules out a proof that assumes the uncorrected finite paired delay products already have the desired sign for all \(t>0\). It does not rule out an estimate on a compact time interval away from zero, or an approximation that treats the remaining factors collectively and preserves the optical gap. It does not settle \(h_T'<h_X'\) in the full parameter family or bound the number of physical capture reversals there.
