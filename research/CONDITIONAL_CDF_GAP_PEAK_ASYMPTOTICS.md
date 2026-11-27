# 条件 CDF 差的概率解释与单峰位置的两端渐近

## 1. 条件分布的两个版本

令
\[
T_c=X_c+Z_c,
\]
其中 \(X_c,Z_c\) 独立，CDF 为 \(F_X,F_T\)，密度为 \(f_X,f_Z,f_T\)。定义
\[
h_X(u)=\frac{f_X(u)}{F_X(u)},\qquad
h_T(t)=\frac{f_T(t)}{F_T(t)},\qquad
q_t(z)=h_X(t-z).
\]
在 \(0<z<t\) 上定义
\[
p_t^-(z)=\frac{F_X(t-z)f_Z(z)}{F_T(t)}
\]
和
\[
p_t^0(z)=\frac{f_X(t-z)f_Z(z)}{f_T(t)}.
\]
它们分别是
\[
\pi_t^-=\mathcal L(Z_c\mid T_c\le t),
\qquad
\pi_t^0=\mathcal L(Z_c\mid T_c=t)
\]
的条件密度。

因为 \(h_T=f_T/F_T\)，有精确 Radon–Nikodym 关系
\[
\boxed{
\frac{p_t^0(z)}{p_t^-(z)}
=\frac{q_t(z)}{h_T(t)}.
}
\tag{1.1}
\]
而
\[
h_T(t)=\int_0^t q_t(z)p_t^-(z)\,dz.
\]
由于 \(h_X\) 严格递减，\(q_t(z)=h_X(t-z)\) 严格递增。因此 \(\pi_t^0\) 相对于 \(\pi_t^-\) 具有严格递增的似然比，特别地
\[
\pi_t^-\le_{\rm lr}\pi_t^0
\]
并且 \(\pi_t^0\) 在一阶随机序意义下向较大的延迟值移动。

## 2. 单峰权重是两个条件 CDF 的精确差

沿用
\[
K_t(x,y)=P_t(x\wedge y)\,[1-P_t(x\vee y)]
\]
及
\[
W_t(y)=\int_0^t A_t(x)K_t(x,y)\,dx,
\qquad A_t(z)=q_t'(z)=-h_X'(t-z)>0.
\]
此前得到
\[
W_t'(y)=p_t^-(y)[h_T(t)-q_t(y)],
\qquad W_t(0+)=W_t(t-)=0.
\]
另一方面，由 (1.1)，
\[
\begin{aligned}
h_T\bigl(P_t^-(y)-P_t^0(y)\bigr)
&=\int_0^y[h_T-q_t(z)]p_t^-(z)\,dz\\
&=W_t(y).
\end{aligned}
\]
因此
\[
\boxed{
W_t(y)=h_T(t)
\left[
\Pr(Z_c\le y\mid T_c\le t)
-\Pr(Z_c\le y\mid T_c=t)
\right].
}
\tag{2.1}
\]
这不仅说明 \(W_t>0\)，也解释了其单峰性：两条件密度的似然比 \(q_t/h_T\) 严格递增，且在唯一点 \(y_t\) 穿过 1。该交点满足
\[
\boxed{h_X(t-y_t)=h_T(t).}
\tag{2.2}
\]
它同时是 CDF 差 (2.1) 的唯一最大点。

## 3. 方差目标变成两个条件分布下的 score 均值差

记
\[
\ell_Z(z)=(\log f_Z(z))'.
\]
此前的精确恒等式为
\[
h_T'-h_X'
=A_t(0)+\operatorname{Cov}_{\pi_t^-}(q_t,\ell_Z).
\]
由 (1.1)，
\[
\operatorname{Cov}_{\pi_t^-}(q_t,\ell_Z)
=h_T\left(
\mathbb E_{\pi_t^0}\ell_Z-
\mathbb E_{\pi_t^-}\ell_Z
\right).
\]
因此
\[
\boxed{
h_T'-h_X'
=A_t(0)-h_T(t)\,D_t,
\qquad
D_t:=
\mathbb E_{\pi_t^-}\ell_Z-
\mathbb E_{\pi_t^0}\ell_Z.
}
\tag{3.1}
\]
所需相对曲率条件等价于
\[
\boxed{
D_t>\frac{A_t(0)}{h_T(t)}.
}
\tag{3.2}
\]
这里 \(D_t\) 是将延迟变量的条件从“总时间不超过 \(t\)”改为“总时间恰等于 \(t\)”时，score 均值的下降量。

还可以直接计算这两个期望：
\[
\mathbb E_{\pi_t^-}\ell_Z=h_T(t),
\qquad
\mathbb E_{\pi_t^0}\ell_Z=\frac{f_T'(t)}{f_T(t)}.
\]
故
\[
\boxed{
D_t=-\frac{h_T'(t)}{h_T(t)}
=-(\log h_T(t))'.
}
\tag{3.3}
\]
于是 (3.2) 与 \(h_T'<h_X'\) 完全等价。条件分布表示的价值在于给出了这个导数比较的概率含义，并为峰值位置和两端渐近提供了可用结构。

## 4. 峰值位置的小时渐近：由光学长度分割

设高频 Robin 反演给出带导数的鞍点渐近
\[
h_X(u)\sim\frac{d_X^2}{4u^2},
\qquad u\downarrow0,
\]
\[
h_T(t)\sim\frac{d_T^2}{4t^2},
\qquad t\downarrow0,
\]
其中 \(d_T=d_X+d_Z\)，\(d_X,d_Z>0\) 分别是 \(X_c,Z_c\) 的高频光学长度。

将其代入峰值方程 (2.2)：
\[
\frac{d_X^2}{4(t-y_t)^2}
\sim
\frac{d_T^2}{4t^2}.
\]
因此
\[
\boxed{
\frac{y_t}{t}\longrightarrow
\frac{d_Z}{d_X+d_Z},
\qquad t\downarrow0.
}
\tag{4.1}
\]
这正是小时间卷积鞍点的光学长度分割点：总光学距离 \(d_T\) 被分成 \(d_X\) 与 \(d_Z\) 两段。

对于当前三层基准，
\[
d_X=\ell_L=\frac7{22},\qquad
d_Z=\ell_R-\ell_L=\frac2{11},\qquad
d_T=\ell_R=\frac12.
\]
故
\[
\boxed{y_t/t\longrightarrow4/11.}
\]

同一渐近还给出
\[
A_t(0)=-h_X'(t)\sim\frac{d_X^2}{2t^3},
\qquad
-h_T'(t)\sim\frac{d_T^2}{2t^3}.
\]
利用 (3.3)，
\[
D_t\sim\frac2t,\qquad
\frac{A_t(0)}{h_T(t)}
\sim\frac{2d_X^2}{d_T^2t}.
\]
所以
\[
\boxed{
D_t-\frac{A_t(0)}{h_T(t)}
\sim
\frac2t\left(1-\frac{d_X^2}{d_T^2}\right)>0.
}
\tag{4.2}
\]
因此对每个具有严格光学间隙的固定 Robin 几何，目标方差不等式在充分小时严格成立。

## 5. 峰值位置的晚时渐近：由指数倾斜矩确定

令
\[
\eta=\lambda_3-c>0,
\]
并沿用最慢极点留数
\[
h_X(t)\sim B_Xe^{-\eta t},
\qquad
h_T(t)\sim B_Te^{-\eta t}.
\]
由已有精确留数乘积，
\[
\boxed{
M:=\frac{B_T}{B_X}
=\mathbb E[e^{\eta Z_c}]>1.
}
\]
将晚时展开代入 (2.2)，得到
\[
B_Xe^{-\eta(t-y_t)}
\sim B_Te^{-\eta t}.
\]
因此
\[
\boxed{
y_t\longrightarrow
y_\infty:=\frac{\log M}{\eta},
\qquad t\to\infty.
}
\tag{5.1}
\]
由严格 Jensen 不等式，
\[
\boxed{
y_\infty>\mathbb E[Z_c].
}
\]
故权重峰值在早时按光学长度比例从原点线性移动；晚时则收敛到由 \(Z_c\) 的指数矩确定的有限延迟尺度。

## 6. 两个条件分布的晚时极限

对 \(\pi_t^-\)，有
\[
p_t^-(z)=\frac{F_X(t-z)f_Z(z)}{F_T(t)}
\longrightarrow f_Z(z),
\]
所以
\[
\pi_t^-\Rightarrow\mathcal L(Z_c).
\]

对 \(\pi_t^0\)，利用
\[
f_X(t-z)\sim B_Xe^{-\eta(t-z)},\qquad
f_T(t)\sim B_Te^{-\eta t},
\]
得到
\[
p_t^0(z)\longrightarrow
\frac{e^{\eta z}f_Z(z)}{M}.
\]
即
\[
\pi_t^0\Rightarrow\mathcal L(Z_c^{(\eta)}),
\]
其中 \(Z_c^{(\eta)}\) 是参数为 \(\eta\) 的指数倾斜分布。

令 \(\ell_Z=(\log f_Z)'\)。由端点渐近，
\[
\mathbb E[\ell_Z(Z_c)]
=\int_0^\infty f_Z'(z)\,dz=0,
\]
而
\[
\mathbb E[\ell_Z(Z_c^{(\eta)})]
=\frac1M\int_0^\infty e^{\eta z}f_Z'(z)\,dz
=-\eta.
\]
因此
\[
\boxed{D_t\longrightarrow\eta.}
\]
同时
\[
\frac{A_t(0)}{h_T(t)}
\longrightarrow\frac{\eta}{M}.
\]
所以晚时的严格裕量为
\[
\boxed{
D_t-\frac{A_t(0)}{h_T(t)}
\longrightarrow
\eta\left(1-\frac1M\right)>0.
}
\tag{6.1}
\]
这从条件分布角度重获了已有的晚时谱留数方差结果。

## 7. 当前能严格关闭的时间区域

对每个固定的、具有严格光学间隙和简单最慢极点的 Robin 参数：

- 小时渐近 (4.2) 证明 (3.2) 在 \(t\downarrow0\) 成立；
- 晚时渐近 (6.1) 证明 (3.2) 在 \(t\to\infty\) 成立。

因此存在依赖于该几何参数的
\[
0<t_-(\kappa)<t_+(\kappa)<\infty
\]
使目标方差不等式在
\[
(0,t_-)\cup(t_+,\infty)
\]
成立。尚未关闭的部分被压缩为有限中间时间区间
\[
[t_-(\kappa),t_+(\kappa)].
\]

对临界参数族 \(\kappa\uparrow\kappa_0\)，要使这成为参数统一结论，还需对小时间鞍点余项、晚时谱余项及最慢谱隙建立一致界；本笔记不把逐参数渐近误写为参数统一证明。

## 8. 下一步的精确目标

中间时间只需证明
\[
\boxed{
D_t=
\mathbb E[\ell_Z(Z_c)\mid T_c\le t]
-\mathbb E[\ell_Z(Z_c)\mid T_c=t]
>
\frac{-h_X'(t)}{h_T(t)}.
}
\]
两个条件分布具有严格递增的似然比，且其 CDF 差正好是单峰权重 \(W_t/h_T\)。下一步应利用这一单峰 CDF 差控制 \(\ell_Z\) 晚时回升造成的正贡献，而不再估计完整 Poincaré 算子范数。

本轮只使用已有 Robin 谱和晚时留数，不重新计算谱根或数值证书。
