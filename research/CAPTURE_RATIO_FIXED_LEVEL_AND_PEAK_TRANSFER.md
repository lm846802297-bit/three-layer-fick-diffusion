# Fixed-level criterion and transfer to the normalized physical capture ratio

This note gives two analytical strengthenings of the ratio reduction. The criterion applies conditionally to any geometry meeting the indicated hypotheses; verification over the entire three-layer family is still open. At the already certified benchmark and its open neighborhood it also yields an ordered peak theorem for the actual normalized capture ratio, not merely for the differentiated auxiliary ratio.

## Normalizations

Write pi_i=P_i(infinity)>0, where P_i are the raw cumulative probabilities. Set U_i=P_i/pi_i and G_i=U_i+U_i'/lambda_1. Let r=G_R'/G_L', Q=G_R/G_L, p=U_R/U_L, gamma=pi_L/pi_R. Thus P_R-P_L=pi_R U_L(p-gamma).

The half-trace spectral formula r=(b/a) M_Z [f_T-(c-a)F_T]/[f_X+(b-c)F_X] is for this normalized r. For the raw filtered derivative ratio, multiply that formula by pi_R/pi_L. They agree at the tie. The previous stability note used raw filtered responses; its peak and zero persistence are unaffected because the multiplicative correction is positive, depends only on parameters, and is 1 at the seed. The comparison with level 1 for raw responses is exactly the comparison with level gamma for normalized responses.

Assume G_L'>0 for all t>0, both normalized filtered responses and both normalized cumulative probabilities tend to 1, and all four vanish at zero. The spectral setting a<c<b with positive remaining half stages guarantees G_L'>0. We use the short-time optical ordering r(0+)=0 and eventual r<0, as occurs in the strict first-inversion setting with longer right optical length. Time functions are smooth for t>0.

## 1. Only a fixed level of r is needed

Let H=G_R-gamma G_L and V=exp(lambda_1 t)(U_R-gamma U_L). Then

    H' = G_L' (r-gamma),
    V' = lambda_1 exp(lambda_1 t) H,
    H(0)=V(0)=0,  H(infinity)=1-gamma.

Suppose {t>0:r(t)>gamma} is empty or a single bounded interval (a_0,b_0). No assumption on r at other levels is made. In these diffusion applications the time functions are analytic; they cannot be equal to the threshold on an interval because their early and late signs are strict. Hence H is strictly decreasing before a_0, increasing on (a_0,b_0), and decreasing after b_0. This statement refers to strict monotonicity and tolerates isolated stationary points at the threshold.

For gamma>1, define h_peak=H(b_0) if the interval exists. If the interval is empty or h_peak<=0, V stays negative for t>0, so there is no reversal. Otherwise H changes from negative to positive once and from positive to negative once. Denote the latter zero by u_2. Then V decreases, increases, and decreases, and tends to -infinity because U_R-gamma U_L tends to 1-gamma<0. It follows that:
- V(u_2)>0: exactly two strict, simple reversals;
- V(u_2)=0: one tangency, no reversal;
- V(u_2)<0: no positive-time zero.

Simplicity of the tangency follows if H'(u_2)<0; in the strict unimodal r application this holds since u_2>b_0 and r<gamma there. If a fixed-level criterion is used in a setting allowing extra isolated touches r=gamma outside (a_0,b_0), then a tangency at u_2 still has no sign change but may have order higher than two. This is a zero-count test at one physically selected level, not a request to prove global log-concavity or a unique peak of r.

For gamma<1, H tends to a positive limit and is eventually decreasing, hence positive on its final branch. The initial negative branch forces the interval to be nonempty, and H changes sign once, from negative to positive. Then V likewise changes sign once, since V tends to +infinity. For gamma=1 the same sign pattern of H follows from its zero limit and strict late decrease. If A_std=lim exp(lambda_1 t)(U_R-U_L) is the signed ground contrast, then V crosses exactly once when A_std>0 and never when A_std<=0. This gives a reusable analytic sufficient condition for all three level cases.

## 2. The peak passes twice, with strict ordering

When r has a single positive component and strictly increases then strictly decreases there, Q is its weighted running average with weight G_L'. Therefore

    Q'= (G_L'/G_L)(r-Q),   [G_L(r-Q)]'=G_L r'.

At the r-maximum, r-Q>0, since strictly smaller historical values have positive weight. At the end of the positive component r=0<Q. On the descending positive portion G_L(r-Q) strictly decreases, hence there is exactly one zero of Q', strictly after the r-peak. Where r stays nonpositive the derivative Q' is strictly negative. If r'<0 on its descending positive portion, then Q''<0 at its peak.

A second weighted average is exact, because U_i solves U_i'+lambda_1 U_i=lambda_1 G_i, with zero initial data:

    p(t)= [integral_0^t exp(lambda_1 s) G_L(s) Q(s) ds]
          /[integral_0^t exp(lambda_1 s) G_L(s) ds].

Let I(t) denote the denominator. Then p'=(I'/I)(Q-p) and [I(Q-p)]'=I Q'. Before and at the Q-peak we have Q>p. After it, I(Q-p) is strictly decreasing, so p has at most one peak. It has a finite peak if it approaches its limit 1 from above, as at the kappa=4 benchmark and its open neighborhood by the positive normalized ground contrast and uniform spectral tail. At its peak p''=(I'/I) Q'<0. In particular, on this neighborhood,

    t_r < t_Q < t_p,       max(r) > max(Q) > max(p) > 1.

The strict height inequalities follow because each peak value of a running average is strictly less than the maximum of its input: the weight is positive on all earlier times. Nondegenerate maxima and the strict derivative signs needed here hold for the fixed-case certificate and its analytical continuation.

For any fixed normalized response pair satisfying this peak theorem, treating gamma as a level gives a sharp classification: 0<gamma<=1 yields one crossing, 1<gamma<max(p) yields two, gamma=max(p) gives a quadratic tangency, and gamma>max(p) gives none. Here p(0+)=0 (from the short optical ordering or the two running averages) and p(infinity)=1 from above. Close to the fold the crossing times obey

    t_± = t_p ± sqrt(2(max(p)-gamma)/(-p''(t_p))) + O(max(p)-gamma).

This last statement is a classification of levels for a fixed response shape. Physically changing the wall parameters changes both gamma and p. A fold in a physical parameter therefore requires a transverse crossing of gamma=max(p) along that path; no such new global parameter fold is asserted here.
