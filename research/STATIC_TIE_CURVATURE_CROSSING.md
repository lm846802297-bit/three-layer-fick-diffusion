# Independent closure of the fixed-benchmark curvature crossing

This note uses the geometry and normalization in STATIC_TIE_CURVATURE_CONTOUR.md. For the untilted suffix density f, let C=(log f)''. Exponential tilting leaves C unchanged. The prior contour certificate proves C<0 on [0.0015,0.010]. The two new scripts establish:
- C'>937000 throughout [0.010,0.012], and the unique zero therein is between 0.01102315539 and 0.01102315546.
- C>0 for all t>=0.012, with a strict positive normalized double-residue margin.

Consequently C has exactly one simple upward zero on [0.0015,infinity). Taken together with the existing short-time estimate C<0 on (0,0.0015], this gives the corresponding all-positive-time result at this benchmark. This package does not independently reprove the short-time estimate.

## Crossing

Use the previous pole-free contour formula for f^(j), now for j=0,1,2,3. The only changes to the quadrature proof are the power seven in both Gaussian bounds and the endpoints 0.010 and 0.012 in the integrable majorant and the time Taylor remainder. All four jets are enclosed on two hundred closed intervals of width 1e-5. Evaluating

    C' = f'''/f - 3 f' f''/f^2 + 2 (f'/f)^3

on the enclosed jets gives the derivative margin. For the final root bracket, direct point enclosures use the same alias and quadrature-tail errors, with no time interpolation. The two contours Q=24 and Q=26 give the same root bracket. The left ground root is refined by bisection within its existing first-root bracket.

## Late tail

Let R_j=k_j^2 be the right-half Robin-Dirichlet eigenvalues, defined by k_j=(2j-1)pi+2 atan(4/k_j). With a=L_1, b=R_1, the residue A_j of V at -R_j is

    A_j = (a/b) (R_j-b)/(R_j-a) * Y_L(-R_j)/Y_R'(-R_j), j>=2.

The signs of the residues follow from the signed derivative Y_R'=sin(k_j/2)(1+24/R_j)/(4k_j). The inverse series follows by moving the q contour to a fixed-width rectangle about the imaginary axis, with horizontal edges at Im(q)=+/-2N pi. The numerator is bounded there on fixed real widths, the right trace is bounded away from zero for large N, and the Gaussian makes the horizontal integrals vanish. The two sides count the symmetric poles twice and hence give the stated s-residues. For t>0, the coefficient bound |A_j|<60 k_j ensures absolute locally uniform convergence with any finite number of time derivatives. Indeed |Y_L(-R_j)|<11+8.8/k_j<12, the rational prefactor is less than one, and |1/Y_R'|<=4 sqrt(k_j^2+16).

Thus f f''-(f')^2 is the double sum of A_i A_j (R_j-R_i)^2 exp(-(R_i+R_j)t) over 2<=i<j. Factor out exp(-(R_2+R_3)t); the (2,3) coefficient is strictly positive. For j<=18 subtract the negative pair terms using their absolute upper bounds at t=0.012. For j>18 bound (R_j-R_i)^2 by R_j^2, extend the inner sum over all i, and bound both absolute one-dimensional sums using |A_j|<60 k_j and (2j-1)pi<k_j<2j pi. The ratio of successive majorants for power p is bounded by ((j+1)/j)^p exp(-8j pi^2 t), decreasing in j. The omitted normalized pair sum does not grow after t=0.012 since every actual pair rate exceeds R_2+R_3. The resulting margin divided by the positive (2,3) coefficient exceeds 0.23913.

## The remaining comparison

For the independent tilted variables X and Z and their sum T, the desired inequality concerns the reversed hazards h_W=F_W'/F_W. The curvature result applies to the density of Z, not the CDFs of X and T. It cannot alone prove h_T'<h_X'. A direct next target on a compact time window is to enclose F_W, f_W and f_W' separately by exact transform inversion for W=X,T. The sign to check, with no division or covariance estimates, is

    (f_T' F_T - f_T^2) F_X^2 - (f_X' F_X - f_X^2) F_T^2 < 0.

For either W the three integrands are generated from the same density transform by multiplication by 1/s, 1 and s, respectively. In particular the CDF integrand is not obtained by separately bounding a lower-limit time integral through t=0. The density curvature zero is not a zero of a risk-rate difference or the physical capture difference.
