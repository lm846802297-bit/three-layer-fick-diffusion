# 第一模唯一逆序：左壁最终占优时的全时间无反超判据

模型及记号沿用 FIRST_INVERSION_LATE_PAIR.md。设
\(a=\mu_{L,1}<b=\mu_{R,1}\) 是唯一逆序，\(p=b/a\)，
\(m=\mu_{R,2}\)，\(\lambda=\lambda_1<a\)，
\(0<\rho=R_L/R_R<1\)。记 \(R=r(-\lambda)\) 为去掉静态比 \(\rho\) 的谱比在全区间基态极点处的值。以下判据只需要 \(\lambda,b,m,\rho,R\)，以及确认唯一逆序；不需要逐时验证。

**定理。** 若
\[
 \frac{\rho b-\lambda}{\rho(b-\lambda)}
 \ \le R\ \le\
 \frac{b(m-\lambda)}{\rho m(b-\lambda)},
\]
则对所有 \(t>0\)，\(\Delta=A_L-A_R>0\)。上、下界之间有严格余量时，它给出开参数区，而不是仅在静态打平、基态抵消的临界面成立。

## 推导

以 \(F=A_L\) 及其归一化拉普拉斯变换为起点。令 \(\nu\) 为所有 \(i\ge2\) 的正常谱对产生的独立等待和 \(Z\) 的分布，\(H=\mathbb E e^{\lambda Z}\)，\(S(t)=\mathbb P(Z>t)\)。这给出
\[
 r(s)=\frac{1+s/a}{1+s/b}\widehat\nu(s),\qquad
 H=\frac{b-\lambda}{b-p\lambda}R .
\]
由全区间谱对右半谱逐项不大于的比较，精确分解
\[
 \widehat F(s)=\frac{F(\infty)}{s}\frac{1+s/b}{1+s/\lambda}\widehat\omega(s)
\]
中的 \(\omega\) 是移去首对因子后的非负概率测度。其任意 \([0,t)\) 都有正概率：有限多阶段一同选零，再用尾和期望趋于零。消去公共 \(1+s/b\) 并反演，得到
\[
 \Delta(t)=\frac{F(\infty)\lambda}{b}(\omega * U)(t),\qquad
 U(t)=A+B e^{-\lambda t}+\rho\Gamma(t),
\]
其中
\[
 A=\frac{b(1-\rho)}{\lambda}>0,\quad
 B=\frac{(b-\lambda)(\rho R-1)}{\lambda},\quad
 \Gamma=pS-(b-p\lambda)e^{-\lambda t}\int_t^\infty e^{\lambda u}S(u)\,du .
\]
为核对反演，可先反演 \((1-\rho r)\widehat F\)，将括号化为
\((1-\rho p)e^{-\lambda t}+A(1-e^{-\lambda t})+\rho pS+
 \rho(b-p\lambda)(e^{-\lambda\cdot}*S)\)，再用
\(\int_0^\infty e^{\lambda t}S(t)\,dt=(H-1)/\lambda\)。

给定 \(Z>t\)，条件化到仍在进行的指数等待阶段。该阶段的剩余时间服从原速率、后续阶段独立，且速率至少为 \(m\)，后续指数矩的子积不大于 \(H\)。因此
\(\mathbb E(e^{\lambda(Z-t)}\mid Z>t)\le mH/(m-\lambda)\)。
用 \(S\) 的积分恒等式整理得 \(\Gamma(t)\ge C S(t)\)，其中
\[
 C=\frac{b-(b-\lambda)R\,m/(m-\lambda)}{\lambda}.
\]
余下不需要单调性：
- 若 \(C\ge0\)，定理下界等价于 \(A+B\ge0\)，故 \(U(t)>0\)（\(A>0,t>0\)）。
- 若 \(C<0\) 且 \(B\le0\)，利用 \(S\le1\) 得
  \(U\ge A+B+\rho C=1-\rho R(b-\lambda)/(m-\lambda)>0\)，
  最后一步由 \(\rho R\le1\) 及 \(m>b\)。
- 若 \(C<0\) 且 \(B>0\)，定理上界等价于 \(A+\rho C\ge0\)，因此 \(U(t)\ge B e^{-\lambda t}>0\)。

卷积为正，得到结论。实际检测 \(R\) 无需乘积截断：
\[
 R=Y_L(-\lambda)/[\rho Y_R(-\lambda)].
\]

## 区间算例与边界

用扩散率 \(2500,1/4,1\)，\(\kappa=7\)，慢层宽分别为 \(0.124\) 与 \(0.12\)，
用区间运算核对了首模、右半第二模、唯一逆序和严格裕量。第一组
\(\rho\approx0.994012\)，左门槛约 \(0.997048\)，\(R\approx0.998177<1\)；
且 \(\rho b/a\approx1.178714>1\)，故只用
\(A_R\le\rho(b/a) A_L\) 无法得出这一结果。
第二组验证 \(R>1\) 的一侧。运行 python first_inversion_residual_band.py。

本判据并不覆盖已有有限时间双零点的参数域。它只排除一个明确的开区内的所有负谷，未排除其补集出现两个彼此分开的负谷，也没有证明第一逆序子族或整个单逆序族的全时间上限。
