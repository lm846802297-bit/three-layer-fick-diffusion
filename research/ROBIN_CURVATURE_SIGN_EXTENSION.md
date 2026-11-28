# Robin score curvature: certified negative and positive outer windows

## 1. Scope and benchmark

This certificate concerns the static equal-resistance three-layer Robin benchmark
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
with right-half width \(1/2\), right diffusivity \(1\), and common Robin parameter
\[
\kappa=4.
\]
Write \(c=\lambda_2\) for the second eigenvalue of the full Robin problem, and
\[
\rho_j=R_j-c,\qquad
A_j=-2C_0 k_j\frac{Y_L(k_j)}{Y_R'(k_j)}
\frac{1-R_j/R_1}{1-R_j/L_1},\qquad k_j=\sqrt{R_j}.
\]
The first half-spectrum roots are distinct and satisfy \(L_1<c<R_1\). The last factor in the residue formula is essential; the common-first-root formula is not applicable to this benchmark.

The curvature numerator is
\[
Q_Z(t)=f_Z(t)f_Z''(t)-f_Z'(t)^2
=\sum_{2\leq i<j}C_{ij}e^{-\lambda_{ij}t},
\]
where
\[
C_{ij}=A_iA_j(\rho_j-\rho_i)^2,
\qquad \lambda_{ij}=\rho_i+\rho_j.
\]
Since \(f_Z(t)>0\), the sign of \(Q_Z\) is the sign of \((\log f_Z)''\).

## 2. Strict negativity on \([0.0015,0.010]\)

The finite sum \(Q_{50}\) uses interval-enclosed roots \(R_2,\ldots,R_{50}\), together with \(L_1,R_1,c\). The root brackets have wave-number radius \(10^{-25}\); interval arithmetic uses 70 decimal digits. The root cells are indexed by the strictly increasing right Prüfer phase. The first left root is indexed by the left Prüfer phase, and the global root \(c\) lies between \(L_1\) and \(R_1\).

The analytic residue estimate gives
\[
|A_j|\leq Bj,\qquad
\pi^2j^2<\rho_j<4\pi^2j^2,
\]
with the computed interval upper bound
\[
B<265.424329828643.
\]
The finite roots also pass the direct interval checks \(|A_j|<Bj\) and \(\rho_j>\pi^2j^2\) for every \(2\leq j\leq50\).

For \(m\geq0\), define
\[
U_m(\tau)=B(4\pi^2)^m\sum_{j\geq2}j^{2m+1}e^{-\pi^2j^2\tau},
\]
and let \(U_{m,N}(\tau)\) denote the same sum restricted to \(j>N\). The double-spectrum estimate is
\[
|Q_Z(t)-Q_N(t)|\leq
2\bigl[S_0(t)S_{2,N}(t)+S_2(t)S_{0,N}(t)\bigr],
\]
so for \(t\geq\tau\),
\[
|Q_Z(t)-Q_N(t)|\leq
2\bigl[U_0(\tau)U_{2,N}(\tau)+U_2(\tau)U_{0,N}(\tau)\bigr].
\]
Each Gaussian series is summed to a finite index and its remaining tail is bounded by a geometric series. The ratio of successive summands is decreasing, since
\[
\frac{r_{j+1}}{r_j}
=\left(\frac{j(j+2)}{(j+1)^2}\right)^{2m+1}e^{-2\pi^2\tau}<1.
\]

At \(\tau=0.0015\) and \(N=50\), the resulting certified error is
\[
\boxed{
\sup_{t\geq0.0015}|Q_Z(t)-Q_{50}(t)|<61.579377296.
}
\tag{2.1}
\]

To bound the finite sum over each box \([a,b]\), put
\[
M_{2,50}(a)=\sum_{2\leq i<j\leq50}
|C_{ij}|\lambda_{ij}^2e^{-\lambda_{ij}a}.
\]
For \(t\geq a\), \(|Q_{50}''(t)|\leq M_{2,50}(a)\). Taylor's theorem therefore gives
\[
Q_{50}(t)\leq Q_{50}(a)+Q_{50}'(a)(t-a)
+\frac12M_{2,50}(a)(t-a)^2.
\]
The right-hand side is convex in \(t-a\), so its maximum on \([a,b]\) is bounded by its two endpoint values. Applying this estimate to 300 equal boxes covering \([0.0015,0.010]\) gives
\[
\max_{\text{300 boxes}}\sup Q_{50}
\leq-5\,517\,850.2501393.
\]
Combining with (2.1),
\[
\boxed{
Q_Z(t)<-5\,517\,788.6707<0,
\qquad 0.0015\leq t\leq0.010.
}
\tag{2.2}
\]

## 3. Strict positivity for \(t\geq0.012\)

The first pair has
\[
C_{23}\in[448153309.5622244431724363,
448153309.5622244431724416]>0,
\]
and
\[
\lambda_{23}\in[334.09888362179124340034159,
334.09888362179124340034161].
\]
For all other pairs, \(\lambda_{ij}>\lambda_{23}\). Define
\[
E(t)=\sum_{(i,j)\ne(2,3)}
\frac{|C_{ij}|}{C_{23}}
 e^{-(\lambda_{ij}-\lambda_{23})t}.
\]
Then
\[
Q_Z(t)\geq C_{23}e^{-\lambda_{23}t}[1-E(t)],
\]
and \(E(t)\) is strictly decreasing.

At \(t=0.012\), the finite sum over \(2\leq i<j\leq30\), excluding \((2,3)\), is enclosed by
\[
E_{\leq30}(0.012)
\in[0.7621564753126098391357796,
0.7621564753126098391357976].
\]
For all omitted pairs, the Gaussian double-spectrum bound gives
\[
E_{>30}(0.012)<9.3441\times10^{-41}.
\]
Consequently,
\[
\boxed{E(0.012)<0.762157<1.}
\]
Monotonicity of \(E\) now yields
\[
\boxed{Q_Z(t)>0,\qquad t\geq0.012.}
\tag{3.1}
\]

## 4. Connection to the existing middle-window certificate

The repository already contains the 70-digit interval certificate for \(N=30\) on \([0.010,0.012]\). It proves
\[
Q_Z'(t)>1.0568131809\times10^9,
\]
with
\[
Q_Z(0.01102)<0<Q_Z(0.01103).
\]
Together with (2.2) and (3.1), this proves that the only possible zero outside the short-time interval is the unique simple zero
\[
\boxed{z_*\in(0.01102,0.01103).}
\]
More precisely, the sign pattern is now rigorously established on
\[
[0.0015,0.010]\cup[0.010,0.012]\cup[0.012,\infty).
\]
The interval \(0<t\leq0.0015\) is closed by the separate short-time reflection-path certificate. Together with the Q30 middle-window certificate and the first-pair late-time bound in Section 3, the current \(\kappa=4\) benchmark has a complete full-time curvature sign partition.

## 5. Correction to the high-frequency normalization

Because the first half-spectrum roots are distinct, removing them from the two canonical products contributes a factor \(L_1/R_1\) at high frequency. With \(q=\sqrt{s-c}\), the exact high-frequency transfer ratio has leading form
\[
\Phi_Z(s)=
\frac{11C_0L_1}{2R_1}e^{-\delta_\ell q}
\left[1+\frac{\gamma_1}{\sqrt{s}}+O(s^{-1})\right],
\qquad
\gamma_1=\frac{\delta_\ell c}{2}-\frac{4\kappa}{5}.
\]
Thus the leading coefficient in the earlier aggregate-tail note must be multiplied by \(L_1/R_1\). The coefficient in the inverse-Laplace density asymptotic must receive the same factor. This constant does not affect \((\log f_Z)''\) or the sign conclusions above.

## 6. Reproducibility

- robin_q50_taylor_box_300.py: source for the new negative-window and late-positive certificates.
- robin_q50_taylor_box_300_input.json: exact benchmark and acceptance parameters.
- robin_q50_taylor_box_300_output.txt: native 70-digit interval output.
- robin_q50_taylor_box_300_provenance.json: SHA-256 fingerprints and runtime metadata.
- robin_q30_interval_certificate.py and its native output: previously archived middle-window certificate.

No full-spectrum recomputation was performed. Only the right roots through index 50 and the three low roots \(L_1,R_1,c\) were enclosed.
