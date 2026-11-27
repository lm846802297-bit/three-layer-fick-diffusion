# \(Z_c\) 密度的对数曲率：谱留数给出的全局障碍

## 1. 对象与问题

沿用首模唯一逆序归约，令
\[
c=\lambda_2,\qquad \rho_i=R_i-c,\qquad \sigma_i=L_i-c,
\qquad 0<\rho_i\le \sigma_i.
\]
指数倾斜后的延迟变量满足
\[
Z_c=\sum_{i\ge2}B_iE_i,
\qquad
B_i\sim\operatorname{Bernoulli}\!\left(1-\frac{\rho_i}{\sigma_i}\right),
\qquad E_i\sim\operatorname{Exp}(\rho_i),
\]
各阶段相互独立。因此
\[
\Phi_Z(s):=\mathbb E e^{-sZ_c}
=\prod_{i\ge2}\frac{1+s/\sigma_i}{1+s/\rho_i}.
\]
严格光学长度差保证 \(\prod_{i\ge2}\rho_i/\sigma_i=0\)，所以无限和没有零点原子并具有连续密度 \(f_Z\)。

本笔记研究 \((\log f_Z)''\)，而不预设 \(f_Z\) 全局对数凹。

## 2. 极点留数

在谱根简单、乘积局部一致收敛的条件下，\(\Phi_Z\) 是亚纯函数，极点位于 \(-\rho_j\)。其留数为
\[
\boxed{
A_j=
\rho_j\left(1-\frac{\rho_j}{\sigma_j}\right)
\prod_{\substack{i\ge2\\i\ne j}}
\frac{1-\rho_j/\sigma_i}{1-\rho_j/\rho_i}.
}
\]
由于 Weyl 渐近给出 \(\rho_i,\sigma_i\asymp i^2\)，而正光学长度差给出垂直带内的高频衰减，Laplace 反演路径可向左移动。对任意 \(\eta\in(\rho_3,\rho_4)\)，有
\[
f_Z(t)=A_2e^{-\rho_2t}+A_3e^{-\rho_3t}
+O(e^{-\eta t}),
\]
且余项的前两阶导数具有相同指数阶控制。

显然，只要 \(\sigma_2>\rho_2\)，就有 \(A_2>0\)。若
\[
\boxed{\rho_2<\sigma_2<\rho_3,\qquad \sigma_3>\rho_3,}
\]
则在 \(A_3\) 的乘积中，\(i=2\) 的分子和分母同时为负，其余因子为正，故
\[
\boxed{A_3>0.}
\]

## 3. 晚时对数凸性定理

**命题。** 若 \(\rho_2<\sigma_2<\rho_3\)、\(\sigma_3>\rho_3\)，并且其余谱因子满足 \(\sigma_i\ge\rho_i>\rho_3\)（\(i\ge4\)），则
\[
\boxed{
(\log f_Z(t))''
\sim
\frac{A_3}{A_2}(\rho_3-\rho_2)^2
e^{-(\rho_3-\rho_2)t}>0,
\qquad t\to\infty.
}
\]

**证明。** 对
\[
f_Z=A_2e^{-\rho_2t}+A_3e^{-\rho_3t}
+O(e^{-\eta t})
\]
直接计算
\[
f_Z''f_Z-(f_Z')^2
=A_2A_3(\rho_3-\rho_2)^2e^{-(\rho_2+\rho_3)t}
+O(e^{-(\rho_2+\eta)t}).
\]
再除以 \(f_Z^2\sim A_2^2e^{-2\rho_2t}\)，即得结论。\(\square\)

因此，在上述谱排序下，\(f_Z\) **不可能在整个正半轴上对数凹**。这个结论针对延迟密度本身，不是对最终 \(r\) 单峰性的反例。

## 4. 已有 Robin 基准满足该条件

已有静态打平基准的严格根区间为
\[
L_2\in[175.7546354113,175.7546354122],
\]
\[
R_3\in[262.1914603604,262.1914603611].
\]
因此 \(L_2<R_3\)，等价于 \(\sigma_2<\rho_3\)。同一根区间和相位排序还给出 \(L_3>R_3\)，即 \(\sigma_3>\rho_3\)。此外 \(L_2>R_2\)，故
\[
\rho_2<\sigma_2<\rho_3,\qquad \sigma_3>\rho_3.
\]
这些区间来自现有归档
[两侧精确特征函数封闭正常谱尾矩](FIRST_INVERSION_STATIC_TIE_UNIFORM_RIGHT_EXACT_MOMENTS.md)，无需重新求谱。

所以，晚时对数凸性已经在实际 Robin 基准中成立，而非仅仅是一般 Bernoulli--指数和的抽象可能性。

## 5. 小时渐近与条件势

令 \(\ell_L,\ell_R\) 分别为左右半区光学长度。由两侧特征函数的高频传递矩阵渐近，
\[
\Phi_Z(s)=C_0
\exp[-(\ell_R-\ell_L)\sqrt{s}]
\left(1+O(s^{-1/2})\right),
\qquad s\to+\infty,
\]
其中 \(C_0>0\)。相应的鞍点反演给出领先小时间项
\[
f_Z(t)\asymp
t^{-3/2}
\exp\left[-\frac{(\ell_R-\ell_L)^2}{4t}\right].
\]
故领先对数曲率为
\[
(\log f_Z(t))''
\sim-\frac{(\ell_R-\ell_L)^2}{2t^3},
\qquad t\downarrow0.
\]
这里的结论是高频鞍点渐近；若要将其用于显式时间区间，仍需补上带导数的统一余项界。

对条件分布 \(\pi_t\)，其负对数密度（忽略归一化常数）为
\[
V_t(z)=-\log F_X(t-z)-\log f_Z(z),
\]
从而有精确恒等式
\[
\boxed{
V_t''(z)=-h_X'(t-z)-(\log f_Z(z))''.
}
\]
因此，Brascamp--Lieb 路线真正需要的不是 \(f_Z\) 全局对数凹，而是相关区域内的比较
\[
\boxed{
-h_X'(t-z)>(\log f_Z(z))''.
}
\]
在 \((\log f_Z)''\le0\) 的区域，因 \(h_X'\le0\)，条件势自动凸；在晚时对数凸区域，则必须由 \(-h_X'\) 提供额外曲率补偿。

若在某个时间区域有 \(B(z)=-(\log f_Z(z))''>0\)，令
\[
A(z)=-h_X'(t-z)\ge0,
\]
则 \(V_t''=A+B\)，且 \(q_t'(z)=A(z)\)，其中 \(q_t(z)=h_X(t-z)\)。一维 Brascamp--Lieb 不等式给出
\[
\operatorname{Var}_{\pi_t}(q_t)
\le\mathbb E_{\pi_t}\frac{A(Z)^2}{A(Z)+B(Z)}.
\]
代入 \(h_T'=-\mathbb E_{\pi_t}A+\operatorname{Var}_{\pi_t}(q_t)\)，得到
\[
\boxed{
h_T'\le
-\mathbb E_{\pi_t}
\frac{A(Z)B(Z)}{A(Z)+B(Z)}<0.
}
\]
这是局部条件方差上界；它本身尚未给出 \(h_T'<h_X'(t)\)，也尚未关闭完整的 \((\log r)''<0\)。

## 6. 文献位置

- Sumita and Masuda, [Classes of probability density functions having Laplace transforms with negative zeros and poles](https://doi.org/10.2307/1427410), *Advances in Applied Probability* (1987)。该文研究具有负实零点和极点的 Laplace 变换所对应的密度类，特别包括指数混合与指数卷积。其分类框架与本问题的谱乘积直接相关；目前可获得的索引文本只有摘要，尚未据此声称其定理已经覆盖本无限谱极限。
- Balabdaoui and Wellner, [Chernoff's density is log-concave](https://doi.org/10.3150/12-BEJ483), *Bernoulli* **20** (2014), 231--244。文中明确使用 \(PF_2\) 与密度对数凹的等价关系，以及 Schoenberg 的 \(PF_\infty\) 判据。当前 \(\Phi_Z\) 的零点和极点顺序不满足直接套用其 \(PF_\infty\) 结构的条件。

## 7. 当前结论与下一步

1. 对 \(Z_c\) 密度的晚时曲率，已得到由前两个正留数控制的严格渐近。
2. 现有静态打平 Robin 基准满足 \(L_2<R_3\)，所以该密度晚时严格对数凸；全局对数凹路线被排除。
3. 高频特征函数给出小时间负曲率的鞍点渐近；显式统一余项尚待完成。
4. 条件势的准确目标是 \(-h_X'(t-z)>(\log f_Z(z))''\)，而不是要求 \(f_Z\) 全局对数凹。

下一步只需研究这个条件曲率补偿在 \(0<t<t_*\) 上是否成立。已有谱数据和数值证书均不需要重算。
