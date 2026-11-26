# 第一半谱模唯一逆序：中间时段的凸下界

沿用 FIRST_INVERSION_RESIDUAL_BAND.md 的模型和记号。假设左壁最终占优
(\(\rho<1\))，并令 \(m=\mu_{R,2},~\ell=\mu_{L,2}>m\)。全区间基态为
\(\lambda=\lambda_1<\mu_{L,1}<\mu_{R,1}=b\)，写 \(R=r(-\lambda)\)。
已有精确反演
\[
 \Delta(t)=\frac{F(\infty)\lambda}{b}(\omega*U)(t),\qquad
 U(t)=A+B e^{-\lambda t}+\rho\Gamma(t),\quad
 A=\frac{b(1-\rho)}\lambda,\quad
 B=\frac{(b-\lambda)(\rho R-1)}\lambda .
\]
这里 \(\omega\) 为正的概率测度，\(\Gamma\) 是由正常半谱对等待之和 \(Z\)
产生的尾项，详见此前笔记。

## 1. 保留当前阶段的零等待概率

正常谱对 \(i\ge2\) 的等待 \(Z_i\) 以概率
\(\mu_{R,i}/\mu_{L,i}\) 为零，否则服从速率 \(\mu_{R,i}\)
的指数律。因此
\[
 H_i=\mathbb E e^{\lambda Z_i}
 =\frac{1-\lambda/\mu_{L,i}}{1-\lambda/\mu_{R,i}},
 \qquad H=\prod_{i\ge2}H_i=\frac{(b-\lambda)R}{b-(b/a)\lambda}.
\]
在 \(Z>t\) 上条件化到按序执行时的活动阶段 \(i\)。剩余指数的矩为
\(\mu_{R,i}/(\mu_{R,i}-\lambda)\)；之后阶段的矩乘积至多 \(H/H_i\)，
而不是 \(H\)。故
\[
 \mathbb E[e^{\lambda(Z-t)}\mid Z>t]
 \le H\sup_{i\ge2}\frac{\mu_{L,i}}{\mu_{L,i}-\lambda}
 = H\frac{\ell}{\ell-\lambda}.
\]
这也适用于无穷级数，因为在 \(Z>t\) 事件上，总有某个有限阶段已越过 \(t\)。
因此旧估计中的 \(m\) 可以换为 \(\ell\)：
\[
 \Gamma(t)\ge C_\ell \mathbb P(Z>t),\qquad
 C_\ell=\frac{b-(b-\lambda)R\ell/(\ell-\lambda)}{\lambda}.
\]

## 2. 保留第二模本身的尾；控制真实曲线而非未卷积项

当 \(C_\ell>0\) 时，因为 \(Z\ge Z_2\)，
\[
 U(t)\ge P(e^{-\lambda t}),\qquad
 P(x)=A+Bx+E x^\gamma,\quad
 E=\rho C_\ell(1-m/\ell)>0,\quad \gamma=m/\lambda>1 .
\]
\(P\) 严格凸。令 \(Y\sim\omega,\ F_\omega(t)=\mathbb P(Y\le t)\) 及
\(x_t=\mathbb E[e^{-\lambda(t-Y)}\mid Y\le t]\)。正卷积加 Jensen 不等式给
\[
 \frac{b\,\Delta(t)}{F(\infty)\lambda F_\omega(t)}\ge P(x_t).
\]
这里 \(F_\omega(t)>0\) 对任意正时间成立。\(x_t\) 连续且不增，从 1 趋到 0。
证明：各 \(\omega\) 的因子可写为 \((E_i-c_i)_+\)，其中 \(E_i\) 为
独立指数变量；有限和的次水平集在 \((E_i,t)\) 中凸，联合指数密度对数凹，
故边缘积分 \(F_\omega\) 对数凹。对无穷等待和取逐点极限保留此性质。
对 \(t>0\)，积分分部给
\[
 x_t=1-\lambda\int_0^\infty e^{-\lambda v}
               F_\omega(t-v)/F_\omega(t)\,dv
\]
（负时间取零）。对数凹使每个固定平移比值不随 \(t\) 下降，故 \(x_t\)
不增。正时间的连续性源自正且对数凹的 CDF；两端极限由
\(e^{-\lambda t}\le x_t\le1\) 与 \(Y<\infty\) 几乎必然成立给出。

因此：若 \(P\) 在 \([0,1]\) 的最小值严格为正，全时间无反超；
若 \(P\) 取负值，则任何严格负谷都须落在同一个候选时间区间
\(\{t:P(x_t)<0\}\) 中。**后一条不证明该区间内只有一个负谷。**

## 3. 区间认证的新参数点

取三层扩散率 \(2500,\ 1/4,\ 1\)，厚度 \(0.38,\ 0.12,\ 0.50\)，
共同 Robin 参数 \(6.2\)。认证只有第一半谱模逆序。
\(R\approx0.97038135\)，低于旧判据左门槛 \(0.98524006\)，
但区间计算给出 \(\min_{[0,1]}P>0.0073\)，故左壁对任何正时间均累计领先。
判据有严格余量，也覆盖该点附近一个开参数区，绝非逐时间采样。
\(B<0\) 时求最小值只需计算
\(x_*=\min\{1,(-B/(E\gamma))^{1/(\gamma-1)}\}\)；否则最小值在 0。
复算：python first_inversion_valley_envelope.py。

旧研究的有限时间相切点在同一固定几何下约为 \(\kappa=5.09022821\)，
仍不满足这个无反超判据。整个第一模唯一逆序区的反超上限仍未证明。
