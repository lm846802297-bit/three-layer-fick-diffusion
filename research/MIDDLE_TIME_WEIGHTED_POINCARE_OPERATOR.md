# 中间时间区域的加权 Poincaré 常数：谱算子、Muckenhoupt 界与 q-适配权重

## 1. 目标与归约

沿用
\[
p_t(z)=\frac{F_X(t-z)f_Z(z)}{F_T(t)},\qquad 0<z<t,
\]
\[
q_t(z)=h_X(t-z),\qquad A_t(z)=q_t'(z)=-h_X'(t-z)>0.
\]
精确条件矩恒等式为
\[
h_T'-h_X'=A_t(0)-\mathbb E_{\pi_t}A_t
+\operatorname{Var}_{\pi_t}(q_t).
\]
所以相对 CDF 曲率 \(h_T'-h_X'<0\) 的目标是
\[
\boxed{
\operatorname{Var}_{\pi_t}(q_t)
<
\mathbb E_{\pi_t}A_t-A_t(0).
}
\tag{W}
\]

定义带权 Poincaré 常数
\[
\boxed{
C_t=
\sup_{\operatorname{Var}_{\pi_t}(g)>0}
\frac{\operatorname{Var}_{\pi_t}(g)}
{\mathbb E_{\pi_t}[g'(Z)^2/A_t(Z)]}.
}
\]
令
\[
\delta_t=1-\frac{A_t(0)}{\mathbb E_{\pi_t}A_t}.
\]
由于 \(q_t'=A_t\)，有
\[
\mathbb E_{\pi_t}\left[\frac{q_t'^2}{A_t}\right]
=\mathbb E_{\pi_t}A_t.
\]
因此
\[
\boxed{C_t<\delta_t}
\tag{P}
\]
是 (W) 的充分条件。这里必须有 \(\delta_t>0\)。

## 2. 精确 Sturm–Liouville 谱问题

相应的自伴算子为
\[
\mathcal L_tg=
\frac1{p_t}
\left(\frac{p_t}{A_t}g'\right)'.
\]
在自然零通量边界条件及相应闭二次型定义下，
\[
C_t=\lambda_1(t)^{-1},
\]
其中
\[
\boxed{
\lambda_1(t)=
\inf_{\substack{\mathbb E_{\pi_t}g=0\\g\ne0}}
\frac{\mathbb E_{\pi_t}[g'^2/A_t]}
{\mathbb E_{\pi_t}[g^2]}.
}
\]
所以中间时间的目标等价于
\[
\boxed{
\lambda_1(t)>
\frac{\mathbb E_{\pi_t}A_t}
{\mathbb E_{\pi_t}A_t-A_t(0)}.
}
\]
这个谱问题不要求 \(p_t\) 对数凹。

还有一个协方差核表示。令
\[
P_t(x)=\int_0^x p_t(z)\,dz,\qquad
\bar P_t(x)=1-P_t(x),
\]
\[
K_t(x,y)=P_t(x\wedge y)\bar P_t(x\vee y).
\]
一维协方差恒等式给出
\[
\operatorname{Var}_{\pi_t}(g)
=\int_0^t\!\!\int_0^t
g'(x)g'(y)K_t(x,y)\,dx\,dy.
\]
令
\[
\boxed{
S_t(x,y)=
K_t(x,y)
\sqrt{\frac{A_t(x)A_t(y)}{p_t(x)p_t(y)}}.
}
\]
经变量替换 \(u=g'\sqrt{p_t/A_t}\)，加权 Poincaré 常数就是该正半定对称核算子的算子范数：
\[
\boxed{C_t=\|S_t\|_{L^2(0,t)\to L^2(0,t)}.}
\]
因此问题可以直接作为一个一维正核谱问题处理。

## 3. 不依赖对数凹性的 Muckenhoupt 上界

设 \(m_t\) 是 \(p_t\) 的中位数。对权重 \(w_t=1/A_t\)，定义
\[
B_t^-=
\sup_{0<x<m_t}
P_t(x)\int_x^{m_t}\frac{A_t(z)}{p_t(z)}\,dz,
\]
\[
B_t^+=
\sup_{m_t<x<t}
\bar P_t(x)\int_{m_t}^{x}\frac{A_t(z)}{p_t(z)}\,dz,
\]
\[
B_t=\max(B_t^-,B_t^+).
\]
一维 Muckenhoupt 判据给出
\[
\boxed{
\frac12 B_t\le C_t\le4B_t.
}
\]
所以
\[
\boxed{4B_t<\delta_t}
\tag{M}
\]
是一个完全显式的充分条件。

为了消去归一化常数，定义
\[
g_t(z)=F_X(t-z)f_Z(z),
\]
\[
G_t(x)=\int_0^x g_t(z)\,dz,\qquad
\bar G_t(x)=\int_x^t g_t(z)\,dz.
\]
则
\[
\boxed{
B_t^-=
\sup_{0<x<m_t}
G_t(x)\int_x^{m_t}\frac{A_t(z)}{g_t(z)}\,dz,
}
\]
\[
\boxed{
B_t^+=
\sup_{m_t<x<t}
\bar G_t(x)\int_{m_t}^{x}\frac{A_t(z)}{g_t(z)}\,dz.
}
\]
这两个表达式只涉及 \(F_X,f_Z,h_X'\) 的一维积分。

### 3.1 端点行为与中间区间的隔离

在严格正光学间隙下，Robin 特征函数的高频鞍点反演给出
\[
f_Z(z)\asymp z^{-3/2}
\exp\left[-\frac{d_Z^2}{4z}\right],
\qquad z\downarrow0,
\]
\[
F_X(u)\asymp u^{1/2}
\exp\left[-\frac{d_X^2}{4u}\right],
\qquad u\downarrow0,
\]
其中 \(d_Z=\ell_R-\ell_L>0\)，而本半区谱乘积的 \(d_X\) 为正光学长度（当前归一化下 \(d_X=\ell_L\)）。同时
\[
A_X(u)=-h_X'(u)\asymp u^{-3}.
\]

对任意固定紧时间区间 \(t\in[t_0,t_1]\subset(0,\infty)\)，上述估计一致成立。于是：

- 当 \(x\downarrow0\)，
  \[
  G_t(x)\asymp x^{1/2}e^{-d_Z^2/(4x)},\quad
  \int_x^{m_t}\frac{A_t}{g_t}
  \asymp x^{7/2}e^{d_Z^2/(4x)},
  \]
  因而 Muckenhoupt 左端乘积为 \(O(x^4)\)。

- 当 \(x\uparrow t\)，令 \(u=t-x\)，则
  \[
  \bar G_t(x)\asymp u^{5/2}e^{-d_X^2/(4u)},\quad
  \int_{m_t}^{x}\frac{A_t}{g_t}
  \asymp u^{-3/2}e^{d_X^2/(4u)},
  \]
  因而右端乘积为 \(O(u)\)。

因此 \(B_t^\pm\) 的上确界不会由端点奇异性产生；它们由严格位于 \(0<x<t\) 内部的区域控制。并且在任意紧时间区间上，\(B_t\) 有统一有限上界。这一步建立了中间时间加权 Poincaré 常数的良定性，不需要 \(p_t\) 对数凹。

## 4. Chen–Wang 变分公式给出的谱适配 Schur 界

Germain–Swan 的 Chen–Wang 变分公式允许对任意严格单调的测试函数 \(H\) 构造 Poincaré 常数上界。对任意正函数 \(r\)，取
\[
H'(z)=-A_t(z)r(z)<0.
\]
代入变分公式，得到
\[
\boxed{
C_t\le
\operatorname*{ess\,sup}_{0<x<t}
\frac{1}{r(x)p_t(x)}
\int_0^t K_t(x,y)A_t(y)r(y)\,dy.
}
\tag{CW}
\]

这可以理解为对正核算子 \(S_t\) 的 Schur 检验。它提供了比单纯使用 \(\sup \tau_t A_t\) 更灵活的测试自由度。

特别地，选用谱指数族
\[
r_\theta(z)=e^{-\theta z},
\qquad 0\le\theta<\eta,\qquad
\eta=\lambda_3-c,
\]
得到
\[
\boxed{
\mathfrak C_{\theta}(t)=
\sup_{0<x<t}
\frac{e^{\theta x}}{p_t(x)}
\int_0^t K_t(x,y)A_t(y)e^{-\theta y}\,dy,
\qquad C_t\le\mathfrak C_\theta(t).
}
\]
因此
\[
\boxed{
\inf_{0\le\theta<\eta}\mathfrak C_\theta(t)<\delta_t
}
\tag{CW-\(\eta\)}
\]
是另一个充分条件。这里 \(\eta\) 取自 Robin 全谱的最慢正常模态；指数测试函数不是任意选取，而是与晚时谱衰减尺度相匹配。

\(\theta=0\) 时，(CW) 的右侧就是下节定义的 \(q_t\)-适配 Stein 因子上确界。

## 5. \(q_t\)-适配 Stein 因子：方差的精确表示

定义
\[
\boxed{
U_t(x)=
\frac{1}{p_t(x)}
\int_x^t(q_t(y)-h_T(t))p_t(y)\,dy.
}
\]
由于 \(q_t\) 单调递增，\(U_t(x)\ge0\)。它满足
\[
(p_tU_t)'=-(q_t-h_T)p_t.
\]
并且
\[
\boxed{
\operatorname{Var}_{\pi_t}(q_t)
=\mathbb E_{\pi_t}[A_t(Z)U_t(Z)].
}
\]
证明只需交换积分次序：
\[
\begin{aligned}
\int_0^t A_t(x)\int_x^t(q_t(y)-h_T)p_t(y)\,dy\,dx
&=\int_0^t(q_t(y)-h_T)p_t(y)\\
&\quad\times[q_t(y)-q_t(0)]\,dy\\
&=\operatorname{Var}_{\pi_t}(q_t).
\end{aligned}
\]

因此目标 (W) 等价于
\[
\boxed{
\mathbb E_{\pi_t}[A_t(Z)(1-U_t(Z))]>A_t(0).
}
\tag{U}
\]
此外，由
\[
\frac{p_t'}{p_t}
=-q_t+(\log f_Z)',
\]
得到一阶方程
\[
\boxed{
U_t'=h_T-q_t+
U_t\bigl(q_t-(\log f_Z)'\bigr),
\qquad U_t(0+)=0,\quad U_t(t-)=1.
}
\]
这把条件方差问题转化成了一个一维边值问题。右端点的极限不是零：由 Robin 小时间渐近，令 \\(u=t-z\\)，则 \\(p_t(z)\\) 与 \\(F_X(u)\\) 同阶，且 \\(U_t(z)=1+O(u^2)\\)，故 \\(U_t(t-)=1\\)。左端点则有 \\(U_t(z)\\sim(h_T-h_X(t))P_t(z)/p_t(z)=O(z^2)\\)。

### 5.1 与最优加权 Poincaré 权重的联系

对任意满足正则性条件的单调中心函数 \(g\)，Heredia–Joulin–Roustant 的最优权重构造为
\[
w_g(x)=
-\frac{1}{g'(x)p(x)}
\int_0^x g(y)p(y)\,dy.
\]
取 \(g=q_t-h_T\)，形式上得到
\[
\boxed{w_{q,t}(x)=\frac{U_t(x)}{A_t(x)}.}
\]
并且
\[
(p_tw_{q,t}q_t')'=-(q_t-h_T)p_t.
\]
在端点正则的截断区间上，定理 3.1 说明相应单调函数是第一非平凡模态。对当前完整条件密度，\(p_t\) 在两个有限端点均趋于零，而该文的定理假设有限端点处密度为正，因此不能直接套用其谱结论。这里无条件成立的是上面的积分恒等式与候选权重能量相等；若要把完整区间的 \(w_{q,t}\) 认定为最优权重，还需对截断区间取极限并验证闭二次型和边界项。

另一方面，如果
\[
A_t(x)>A_t(0)
\]
且
\[
\boxed{
U_t(x)\le1-\frac{A_t(0)}{A_t(x)}
}
\tag{U-pointwise}
\]
在整个区间成立，那么
\[
\frac{A_t-A_t(0)}{A_t^2}
\ge\frac{U_t}{A_t}=w_{q,t}.
\]
于是目标权重支配 \(q_t\)-适配权重，直接推出 (W)。

这是一个较强的点态充分条件；若它在某些位置失效，并不意味着平均不等式 (U) 失效。

## 6. 本轮得到的结果与真正剩余的工作

已经建立：

1. 中间时间方差问题等价于一个明确的加权 Poincaré 谱隙问题。
2. Muckenhoupt 判据给出完全由 Robin 谱函数构成的双侧积分上界。
3. 高频端点渐近保证这些积分有限，且中间时间的 Poincaré 常数在任意紧时间区间上统一有界。
4. Chen–Wang 变分公式给出可优化的正核 Schur 上界；谱指数族 \(e^{-\theta z}\) 提供一个与 Robin 慢模相匹配的测试族。
5. \(q_t\)-适配 Stein 因子 \(U_t\) 给出条件方差的精确一维表示。

**尚未证明：**
\[
C_t<\delta_t
\]
在整个所需中间时间区间成立。端点有限性不等于目标谱隙下界；本轮没有把这一点写成已完成的物理定理。

下一步应集中证明下列二者之一：

\[
\inf_{0\le\theta<\eta}\mathfrak C_\theta(t)<\delta_t,
\]
或者直接证明平均条件
\[
\mathbb E_{\pi_t}[A_t(1-U_t)]>A_t(0).
\]

后者更弱、更接近实际需要；前者则给出一个完整的加权 Poincaré 常数控制。

## 7. 可复用的前人理论

- A. Saumard (2019), *Weighted Poincaré inequalities, concentration inequalities and tail bounds related to Stein kernels in dimension one*, Bernoulli 25, 1692–1726. [DOI](https://doi.org/10.3150/19-BEJ1117). 给出任意连通支撑一维密度的 Stein-kernel 加权 Poincaré 不等式。
- G. Germain and Y. Swan (2023), *A note on one-dimensional Poincaré inequalities by Stein-type integration*, Bernoulli 29, 1714–1740. [DOI](https://doi.org/10.3150/22-BEJ1518). 给出 Chen–Wang 变分公式、Stein kernel 上下界及加权 Poincaré 常数的迭代估计。
- D. Heredia, A. Joulin and O. Roustant (2025), *On one dimensional weighted Poincaré inequalities for Global Sensitivity Analysis*, Journal of Mathematical Analysis and Applications. [DOI](https://doi.org/10.1016/j.jmaa.2025.129992). Theorem 3.1 根据给定的单调函数构造使其饱和不等式的最优权重。

这些都是已有数学工具。本项目的特定工作是将它们用于 Robin 谱诱导的条件密度 \(p_t\)，并把谱隙阈值与实际反向风险率曲率联系起来。

## 8. 计算状态

本轮是解析推导与文献核对，没有启动新的谱根计算或数值扫描。已有输入、谱根证书和晚时渐近结果均直接沿用。


## 9. 直接协方差恒等式：无需估计完整 Poincaré 谱隙

上一节的目标是为所有测试函数控制 \(C_t\)，但实际只需要控制特定函数 \(q_t\)。利用 Robin 延迟密度的正光学间隙，
\[
f_X(0+)=f_Z(0+)=0,
\]
可得到一个更直接的恒等式。

### 9.1 两个分部积分恒等式

记
\[
\ell_Z(z)=(\log f_Z(z))'.
\]
由于
\[
p_t(z)=\frac{F_X(t-z)f_Z(z)}{F_T(t)},
\qquad
\partial_t\log p_t(z)=q_t(z)-h_T(t),
\]
并且端点项因 \(F_X(0)=f_Z(0)=0\) 消失，有
\[
\boxed{
\mathbb E_{\pi_t}[q_t(Z)]
=
\mathbb E_{\pi_t}[\ell_Z(Z)]
=
h_T(t).
}
\tag{9.1}
\]
第二个等式也可直接由
\[
\int_0^tF_X(t-z)f_Z'(z)\,dz
=
\int_0^tf_X(t-z)f_Z(z)\,dz
=f_T(t)
\]
得到。

对 \(h_T=\mathbb E_{\pi_t}\ell_Z\) 关于 \(t\) 求导，得到
\[
\boxed{
h_T'(t)
=
\operatorname{Cov}_{\pi_t}
\bigl(q_t(Z),\ell_Z(Z)\bigr).
}
\tag{9.2}
\]
所以
\[
\boxed{
h_T'-h_X'
=
A_t(0)+
\operatorname{Cov}_{\pi_t}
\bigl(q_t(Z),\ell_Z(Z)\bigr).
}
\tag{9.3}
\]
因此，目标 (W) 等价于
\[
\boxed{
\operatorname{Cov}_{\pi_t}
\bigl(q_t(Z),(\log f_Z)'(Z)\bigr)
<-A_t(0).
}
\tag{9.4}
\]

这把条件方差比较改写为一个反向风险率与延迟密度 score 之间的协方差比较。

### 9.2 正核的单峰结构

令
\[
K_t(x,y)=P_t(x\wedge y)\bar P_t(x\vee y)
\]
以及
\[
\boxed{
W_t(y)=\int_0^t A_t(x)K_t(x,y)\,dx.
}
\]
Hoeffding–Stein 协方差恒等式给出
\[
\operatorname{Cov}_{\pi_t}(q_t,\ell_Z)
=
\int_0^t\!\!\int_0^t
A_t(x)\ell_Z'(y)K_t(x,y)\,dx\,dy,
\]
从而
\[
\boxed{
h_T'-h_X'
=
A_t(0)+
\int_0^tW_t(y)(\log f_Z(y))''\,dy.
}
\tag{9.5}
\]

这里 \(W_t\) 是严格正的权重。进一步将其拆开：
\[
\begin{aligned}
W_t(y)
={}&\bar P_t(y)\int_0^y A_t(x)P_t(x)\,dx\\
&+P_t(y)\int_y^t A_t(x)\bar P_t(x)\,dx.
\end{aligned}
\]
直接求导并使用
\[
\int_0^t A_t(x)\bar P_t(x)\,dx
=
\mathbb E_{\pi_t}[q_t]-q_t(0)
\]
得到一个更强的精确公式：
\[
\boxed{
W_t'(y)=p_t(y)\,[h_T(t)-q_t(y)].
}
\tag{9.6}
\]

由于 \(q_t\) 严格递增，且
\[
q_t(0)<h_T(t)<q_t(t-)=+\infty,
\]
存在唯一 \(y_t\in(0,t)\) 满足
\[
\boxed{q_t(y_t)=h_T(t).}
\]
因此
\[
W_t'(y)>0\quad(0<y<y_t),\qquad
W_t'(y)<0\quad(y_t<y<t).
\]
也就是说，**条件方差所诱导的 score-curvature 权重 \(W_t\) 严格单峰**。其峰值位置不是任意拟合参数，而由反向风险率方程精确确定。

结合分部积分，
\[
\operatorname{Cov}_{\pi_t}(q_t,\ell_Z)
=
-\int_0^tW_t'(y)\ell_Z(y)\,dy.
\]
又因 \(\int_0^tW_t'=0\)，可减去 \(\ell_Z(y_t)\)：
\[
\boxed{
\operatorname{Cov}_{\pi_t}(q_t,\ell_Z)
=
-\int_0^tW_t'(y)
[\ell_Z(y)-\ell_Z(y_t)]\,dy.
}
\tag{9.7}
\]

这给出了一个直接的符号判据：如果 score \(\ell_Z\) 在 \(y_t\) 左侧严格高于其在 \(y_t\) 的值、右侧严格低于该值，那么协方差严格为负。当前 Robin 密度的 score 在晚时会因对数凸尾部而重新上升，故还需要控制其越过 \(\ell_Z(y_t)\) 后的尾部正贡献。

### 9.3 与 Stein 因子的关系及端点

由协方差核表示，
\[
W_t(y)=p_t(y)U_t(y).
\]
因此 \(W_t\) 单峰等价于
\[
(p_tU_t)'=p_t(h_T-q_t).
\]
端点小时间渐近给出
\[
U_t(z)=O(z^2),\qquad z\downarrow0,
\]
以及
\[
U_t(z)=1+O((t-z)^2),\qquad z\uparrow t.
\]
所以
\[
W_t(0+)=W_t(t-)=0,
\]
且 (9.5) 中的加权曲率积分收敛。

### 9.4 新的中间时间目标

由 (9.5)，真正需要证明的条件方差补偿不再要求估计完整算子范数 \(C_t\)，而是只需证明
\[
\boxed{
\int_0^tW_t(y)(\log f_Z(y))''\,dy<-A_t(0).
}
\tag{9.8}
\]
这比 \(C_t<\delta_t\) 更贴近所需结论，因为它只针对实际测试函数 \(q_t\)，而不是对所有 Sobolev 函数作统一控制。

已知
\[
(\log f_Z)''<0
\]
在小时区间成立，而特定 Robin 谱排序下它在晚时转正。式 (9.8) 将中间时间问题精确化为：由 \(q_t(y_t)=h_T(t)\) 定位的单峰权重，能否使早时负曲率的积分贡献严格压过晚时正曲率及边界项 \(A_t(0)\)。

**本轮已完成精确恒等式和权重单峰性；式 (9.8) 的全中间时间严格估计仍待证明。** 这一步不需要重算晚时留数，也不需要重新证明 \(Z_c\) 的全局对数凹性。
