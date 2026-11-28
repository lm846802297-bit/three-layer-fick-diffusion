# Pole-free contour certificate for the static-tie curvature gap

This is an independent bridge for the fixed three-layer benchmark, not a theorem for the full parameter family. The layer widths are 25/66, 4/33, 1/2; diffusivities are 25, 1/4, 1; and the Robin coefficient is 4. Write V for the Laplace transform of the suffix distribution after both first half-spectrum factors have been removed. Its density f satisfies

    f(t) > 12.70,     (log f)''(t) < -2100,     0.0015 <= t <= 0.010.

The density lower bound pertains to the untilted, normalized suffix. The curvature bound also holds after any admissible normalized exponential tilt, since its log density differs by an affine function of time. The only spectral numbers computed are the two first half-roots required by the definition of V. No higher eigenvalues, residues, finite spectral curvature polynomial or previously certified short-time remainder enter this bridge.

## 1. Exact transfer functions and inversion

Put q = sqrt(s), l = 5/66, n = 8/33, v = 1/5, u = 2, kappa = 4, and define

    Y_L(q^2) = cosh(n q) cosh(l q) + 10 sinh(n q) sinh(l q)
             + (kappa/q)[ v cosh(n q) sinh(l q) + u sinh(n q) cosh(l q) ],
    Y_R(q^2) = cosh(q/2) + (kappa/q) sinh(q/2),
    V(q^2) = (Y_L/Y_R)(q^2) (1+q^2/b)/(1+q^2/a).

These expressions follow by multiplying the two propagation matrices on the left for state (y,D y'/q), starting at (1,kappa/q); the right half needs only one matrix. At s=0 the two traces agree, by equality of the left and right resistance integrals. The first left wave number is bracketed in the script by the increasing, continuously lifted phase

    atan2(10 sin(d), cos(d)) + n k,   d = l k - atan2(kappa v,k),

at phase pi/2 (there is no branch change up to the bracket). The first right wave number solves k = pi + 2 atan2(kappa,k); its scalar residual is strictly increasing. No other roots are needed.

For q=x+iy with x>=12, the transfer-bound factors are

    |Y_L| <= (11+8.8/x) exp((l+n)x),
    |Y_R| >= [(1-4/x)-(1+4/x)exp(-x)] exp(x/2)/2,
    |(1+q^2/b)/(1+q^2/a)| <= (a/b)[1+(b-a)/x^2].

The last bound uses |q +/- i sqrt(a)| >= x. Thus |V(q^2)| < 40 exp(-2x/11), and there are no singularities in the right-hand region used. Upon deforming the Bromwich line to q=Q+iy,

    f^(j)(t) = integral_R G_j(y,t) dy,
    G_j(y,t) = (Q+iy)^(2j+1) V((Q+iy)^2) exp(t (Q+iy)^2)/pi,
    j=0,1,2.

For the deformation, use the Bromwich line Re(s)=Q^2. In q coordinates its ends can be connected horizontally at Im(q)=+/-Y from Q to sqrt(Q^2+Y^2). On the part with Re(q)<=Y/2, exp(t Re(q^2)) gives Gaussian decay, and on the rest exp(-2 Re(q)/11) wins over the polynomial prefactors because Re(q^2)<=Q^2. This also justifies differentiating the deformed integral with respect to t.

## 2. Two explicit quadrature errors

Apply the unit-step trapezoid rule to G_j. It is analytic for |Im(y)|<=12. On each strip boundary the integral of its absolute value is bounded by

    U = (80/pi) exp(0.010 (Q+12)^2)
        integral_0^infty (Q+12+y)^5 exp(-0.0015 y^2) dy.

Here |q|>=12 so power 5 dominates all three derivative weights; for this estimate we discard the factor exp(-2 Re(q)/11). The moments in U obey the elementary Gaussian recurrence used in the script. Periodize G_j on a unit interval and expand in Fourier series. For each nonzero Fourier integer m, shift the contour by 12 in the decaying direction. The coefficient is then bounded by U exp(-24 pi |m|). Absolute summability justifies evaluation of the Fourier series at zero, giving alias error at most 2 U r/(1-r), r=exp(-24 pi).

Keep only sample indices -260 through 260. For n>=261 bound each summand by (40/pi) exp(0.010 Q^2) (Q+n)^5 exp(-0.0015 n^2). The ratios in this majorant are bounded by (262/261)^5 exp(-0.0015*523) < 1. Sum the geometric bound on both sides. This is a tail of contour quadrature, not a tail of a spectral series.

## 3. Full time intervals, not samples

On each of 850 cells of width 1e-5, expand the finite trapezoid sum through degree four in time about the midpoint. For its j-th derivative each term has the form z p^j exp(pt) with p=(Q+iy)^2. The time Taylor remainder is bounded by

    |z| |p|^(j+5) max(exp(0.0015 Re p),exp(0.010 Re p)) (5e-6)^5/5! .

The script sums this remainder, adds both quadrature errors, and checks an interval enclosure of f''/f - (f'/f)^2 on the entire cell. The two contours Q=24 and Q=26 independently pass these checks. The samples and all interval constants are computed from the exact rational geometry (or outward-enclosed decimal grid endpoints); no fitted constants enter the test. The proof is computer-assisted: the analytic quadrature and interpolation bounds are explicit and the finite checks use 55-decimal-digit interval arithmetic.

Run: python static_tie_curvature_contour.py 24; python static_tie_curvature_contour.py 26. Only mpmath is required.

## Scope

This closes the stated sign gap for the fixed benchmark. If the other three windows are taken from the same benchmark and their certificates are verified, it supplies the missing left-side bridge for the one-zero curvature statement. It does not itself verify those other windows, establish a hazard-derivative comparison, prove unimodality of the original ratio, or prove a full-family reversal bound.
