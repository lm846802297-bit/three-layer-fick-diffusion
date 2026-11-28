# Direct CDF-derivative comparison on the independent wall box

Use the fixed layers and release point of REVERSED_HAZARD_FIRST_FACTOR_TAIL_BOX.md, with independent \(\kappa_L,\kappa_R\in[3.9,4.1]\) and \(\eta=\lambda_3-\lambda_2\). The certificate now gives
\[
 e^{\eta t}(h_X'-h_T')>66,\qquad t\ge0.062.
\]
The prior short-time result on \(0<t\le0.0015\) is unchanged; the interval between those two pieces is not covered by this argument. The small-box all-time theorem also remains unchanged.

## Reduce the nonlinear comparison before bounding tails

Use normalized CDFs, so \(0<F_T\le F_X\le1\) by the independent nonnegative delay \(T=X+Z\). The ordered suffix and log-concavity of \(F_X\) also give \(h_T\ge h_X\ge0\): \(h_T(t)\) is a conditional average of \(h_X(t-Z)\), and \(h_X\) is nonincreasing. The spectral ordering proof in the preceding wall-box note applies on this box; the phase and low-mode bracket checks needed for it are included here.

At any time for which \(F_T''<0\) and \(F_X''-F_T''>0\),
\[
\begin{aligned}
h_X'-h_T'
&=\frac{F_X''}{F_X}-\frac{F_T''}{F_T}+h_T^2-h_X^2\\
&\ge\frac{F_X''-F_T''}{F_X}
\ge F_X''-F_T''>0.
\end{aligned}
\]
The sufficient conditions concern only two linear observables. In particular, the error cost of normalizing by each CDF and squaring two separately bounded density ratios is no longer needed. No sign condition on \(\psi_X'\) or \(f_Z'/f_Z\) is used.

## The joint transform and two explicit poles

Put \(p=s-c\), \(c=\lambda_2\), and \(q=\sqrt p\). In the usual transfer notation the two unnormalized CDF transforms are \(Y_R(q)(p+\lambda_1)/[qD(q)(p+b)]\) and \(Y_L(q)(p+\lambda_1)/[qD(q)(p+a)]\). The normalization constants \(N_X,N_T\) are their positive residues at \(s=0\). Thus the transform of \(F_X''-F_T''\) is exactly
\[
 \frac{s^2(p+\lambda_1)}{qD(q)}
 \frac{Y_R(q)(p+a)N_T-Y_L(q)(p+b)N_X}
      {(p+b)(p+a)N_X N_T}.
\]
There are no initial-value terms because the CDFs and densities are flat at zero. The factor \(s^2\) removes the stationary pole. The script checks both this shared-numerator enclosure and the triangle bound on each central contour cell. On this wide box the triangle estimate is the sharper of the two; using a shared denominator does not by itself eliminate interval dependency.

Shift to the contour with central \(\operatorname{Re}p=-168\) and outer rays at \(\operatorname{Re}q=8\). The split spectrum below 168 has exactly three modes \(L_1,R_1,R_2\): the phase check gives \(L_2>170\), and \(R_3>(5\pi)^2\). Interlacing therefore allows at most four full modes below 168. The existing brackets and the additional fourth-root bracket account for all four. After stripping \(\lambda_1\) and removing the stationary pole, the two retained poles are \(\lambda_3\) and \(\lambda_4\). The negative fourth-mode coefficient in the contrast must actually be included; it cannot be discarded as favorable.

For the leading survival gap \(K_T/K_X-1\), keep factors \(i=2,3\) in the delay moment product at \(s=-\eta\). Use \(L_2>170\) and \(R_2\) bracketed for \(i=2\); for \(i=3\), the optical phase bound gives \(L_3>380\) and the right Robin equation gives \(R_3<270\). The factor excess is increasing in the left eigenvalue and decreasing in the right eigenvalue with the other inputs fixed. This proves a relative gap \(>0.37929\) throughout the box.

## Verification

For each linear observable the retained fourth-mode negative part, after scaling by \(e^{\eta t}\), is bounded from below by its negative coefficient times \(e^{-(\lambda_4-\lambda_3)\tau}\). The contour remainder is bounded by its Gaussian-strip mass times \(e^{-(168-\lambda_3)\tau}\); both bounds hold on the entire half-line \(t\ge\tau\). At \(\tau=.062\), the resulting scaled lower bounds are \(66.7428\ldots\) for \(F_X''-F_T''\) and \(1261.4334\ldots\) for \(-F_T''\). This certifies the claimed sign without any time grid. At \(.060\) this particular lower bound for the contrast is negative; that is an inconclusive bound, not a counterexample and not evidence of a sign change.

Run the companion script with no arguments using mpmath 1.4.1. Only the stated default parameter box and the positive reported margins are claimed.
