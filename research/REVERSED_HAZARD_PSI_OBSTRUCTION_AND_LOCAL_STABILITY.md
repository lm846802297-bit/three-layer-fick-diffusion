# Why the density-derivative shortcut fails, and an open-neighborhood reverse-hazard theorem

Notation: for the tilted suffix distributions X and T, let h_W=(log F_W)', psi_X=f_X'/F_X. When T=X+Z independently with smooth density f_Z, condition on X+Z<=t and put U=t-Z. Write w_t for the conditional CDF of U, and H_t=h_T(t)+h_X(t). This note checks the exact one-variable kernel identity and proves an analytic parameter extension of the existing all-time comparison. It does not prove the comparison for the whole first-inversion family or a fixed-level valley exclusion for the filtered ratio r.

## 1. Exact identity and a rigorous obstruction

The convolution for F_T and two differentiations give E[h_X(U)]=h_T and E[psi_X(U)]=f_T'/F_T. The zero-time flatness of the CDF and density removes boundary terms. Therefore

    h_X'-h_T' = (h_T+h_X)(h_T-h_X) + psi_X(t)-E[psi_X(U)]
              = integral_0^t w_t(u) [psi_X'(u)-H_t h_X'(u)] du.

The boundary terms at u=0 vanish as well: F_X and w_t decay exponentially in 1/u, whereas the ratios below have polynomial growth.

The exact transfer formula and the same saddle contour as in STATIC_TIE_REVERSE_HAZARD_CONTOUR.md yield, for either W,

    F_W(t)=C_W exp(c t-d_W^2/(4t)) t^(-1/2) J_W(t),
    J_W=1+O(t),  J_W^(k)=O(t^(1-k)), k=1,2,3.

Here d_X is the left optical length and d_T is the right optical length. These are differentiated estimates, not a formal termwise differentiation. In fact, on q=d_W/(2t)+iv/sqrt(t), J_W is a fixed Gaussian average of a normalized analytic transfer ratio A(q). On Re(q)=x sufficiently large, uniformly in Im(q), A=1+O(x^-1); Cauchy's estimate gives A^(m)=O(x^(-m-1)). Each t derivative q^(ell) is O(t^(-ell-1)(1+sqrt(t)|v|)), so a chain-rule term with m factors has net order t^(1-k) times a fixed polynomial in |v|. Gaussian integrability justifies the differentiation. All constants are uniform on compact neighborhoods with positive layer lengths and diffusivities, because the leading transfer coefficient stays away from zero and the reflected exponential lengths stay bounded below.

Putting A_X=d_X^2/4, one obtains h_X=A_X/t^2-1/(2t)+c+O(1), h_X'=-2A_X/t^3+1/(2t^2)+O(t^-1), and h_X''=6A_X/t^4-1/t^3+O(t^-2). Since psi_X=h_X'+h_X^2,

    psi_X'(t) = -4 A_X^2/t^5 + O(t^-4).

Thus psi_X'(t)>=0 for every t is impossible, including at the already settled baseline. More strongly, for any fixed t>0 the complete pointwise bracket psi_X'(u)-H_t h_X'(u) is strictly negative for sufficiently small u. Positivity of the integral cannot follow from a globally nonnegative bracket.

Late in time, if eta is the common leading decay rate and K_X>0 is the leading survival amplitude, the density derivative series gives psi_X'(t)=eta^3 K_X exp(-eta t)+o(exp(-eta t))>0. The short- and long-time signs of psi_X' therefore differ.

## 2. A bounded short-slice deduction

Let Phi_t(u)=psi_X(u)-H_t h_X(u). For a cutoff 0<a<t, integrate the short slice by parts using w_t'=F_X(u) f_Z(t-u)/F_T(t). This gives exactly

    integral_0^a w_t Phi_t'
      = w_t(a) Phi_t(a) -
        [integral_0^a (f_X'(u)-H_t f_X(u)) f_Z(t-u) du]/F_T(t).

If f_X' is nonnegative on (0,a) and Phi_t(a)>=0, the short slice is bounded below by

    - f_X(a) sup_{v in [t-a,t]} f_Z(v) / F_T(t).

On any positive compact time window, the cutoff can be chosen uniformly small: f_X'/f_X tends to +infinity at zero, whereas H_t is bounded there in the exterior t variable. The penalty is exponentially small in 1/a, times a polynomial. A sufficient inequality on the remaining slice is that its integral exceed this penalty. This is not a claim that the remaining slice has one sign.

## 3. Conditional endpoint theorem and a genuine parameter extension

For any local family with positive finite layer data, simple separated full eigenvalues, a<c<b<L_2 (a=L_1, b=R_1, c=lambda_2), and d_T>d_X, the saddle formula above gives

    h_X'-h_T' = (d_T^2-d_X^2)/(2t^3)+O(t^-1)>0

for short time, uniformly on a sufficiently small compact parameter neighborhood. The normalized CDF transforms are the same formulas as in the fixed-case note with the perturbed half traces and eigenvalues. They are probability transforms: each suffix pairing has lambda_{i+1}<=R_i and lambda_{i+1}<=L_i for i>=2, by inserting the other half's ground eigenfunction in each min-max trial space. This only needs L_1<R_1<L_2; it does not require the separate ordering R_i<=L_i for the delay Z.

For late time assume the common leading survival terms of X,T are K_X exp(-eta t), K_T exp(-eta t), with eta=lambda_3-c and K_T>K_X>0, and the next-rate gap is positive. Then

    h_X'-h_T' = eta^2 (K_T-K_X) exp(-eta t)
                + O(exp(-(eta+delta)t))

for some delta>0. This holds uniformly on a sufficiently small compact neighborhood with a positive lower bound on the two strict gaps (rate gap and amplitude difference). For example, shift the exact CDF inverse Laplace contour to a fixed real line between the third and fourth full poles in p=s-c. The pole separation persists locally; on the compact portion the characteristic determinant is uniformly separated from zero, and the optical estimate at large imaginary p provides a common integrable envelope even after multiplication by the first two time-derivative weights. Normalization and residues vary continuously by simplicity. If T=X+Z, then K_T/K_X=E exp(eta Z)>1 whenever the delay is nontrivial and has a moment beyond eta; at the benchmark the strict amplitude gap is also independently certified directly by the two residues.

At the fixed kappa=4 benchmark, all required low-mode and amplitude inequalities are strict and the all-time comparison was already certified in STATIC_TIE_REVERSE_HAZARD_CONTOUR.md and STATIC_TIE_REVERSE_HAZARD_TAIL.md. Fix common short and late cutoffs from the two uniform endpoint proofs. Its strictly positive h_X'-h_T' on the intervening compact window has a positive minimum. Joint continuity of the CDF jets in t and all local layer and wall parameters follows from their exact inverse Laplace formulas and the same optical domination on any positive compact t window; both CDFs are positive there. Shrinking the parameter box preserves the minimum sign. Thus the full-time reversed-hazard comparison holds throughout a single open neighborhood in layer thickness, all three diffusivities, and both independently varying wall reaction coefficients (release position remains at the middle interface).

The parameter extension uses no new scan; the compact-window seed still relies in part on the existing interval proof. It gives no explicit neighborhood size and does not make psi_X monotone or establish a global two-component exclusion for r.
