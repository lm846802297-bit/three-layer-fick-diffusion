# Robin 谱更新核与反向风险率导数比较

## 1. 研究目标

沿用
\[
a=L_1<c=\lambda_2<b=R_1,
\qquad
L_i\ge R_i\ge\lambda_{i+1},\quad i\ge2,
\]
以及指数倾斜后的独立随机变量
\[
T_c=X_c+Z_c.
\]
记
\[
h_X(t)=\frac{f_X(t)}{F_X(t)},\qquad
h_T(t)=\frac{f_T(t)}{F_T(t)}.
\]

条件方差目标
\[
\operatorname{Var}_{\pi_t}(h_X(t-Z))
<
\mathbb E_{\pi_t}[-h_X'(t-Z)] + h_X'(t)
\]
等价于
\[
\boxed{h_T'(t)<h_X'(t).}
\tag{1.1}
\]
本笔记将 score-curvature 积分和 Robin 谱乘积分别改写成更直接的风险率与更新核对象。

## 2. 单峰权重积分的精确消去

定义
\[
p_t^-(z)=\frac{F_X(t-z)f_Z(z)}{F_T(t)},\qquad
p_t^0(z)=\frac{f_X(t-z)f_Z(z)}{f_T(t)}.
\]
它们分别是
\[
\mathcal L(Z_c\mid X_c+Z_c\le t),
\qquad
\mathcal L(Z_c\mid X_c+Z_c=t)
\]
的条件密度。令
\[
q_t(z)=h_X(t-z),\qquad \ell_Z(z)=(\log f_Z(z))'.
\]
则
\[
\frac{p_t^0(z)}{p_t^-(z)}
=\frac{q_t(z)}{h_T(t)}.
\tag{2.1}
\]
并且
\[
h_T(t)=\mathbb E_{\pi_t^-}[\ell_Z(Z)].
\]
因为
\[
\partial_t\log p_t^-(z)=q_t(z)-h_T(t),
\]
对上式求导得到
\[
\boxed{
h_T'(t)=
\operatorname{Cov}_{\pi_t^-}(q_t(Z),\ell_Z(Z)).
}
\tag{2.2}
\]

令
\[
W_t(z)=h_T(t)\left[
F_{\pi_t^-}(z)-F_{\pi_t^0}(z)
\right].
\]
由 (2.1)，
\[
W_t'(z)=p_t^-(z)[h_T(t)-q_t(z)],
\qquad W_t(0+)=W_t(t-)=0.
\]
一维协方差核恒等式和分部积分给出
\[
\operatorname{Cov}_{\pi_t^-}(q_t,\ell_Z)
=\int_0^tW_t(z)\ell_Z'(z)\,dz.
\]
因此得到精确恒等式
\[
\boxed{
\int_0^tW_t(z)(\log f_Z(z))''\,dz=h_T'(t).
}
\tag{2.3}
\]

又因 \(h_X'(t)=-A_t(0)\)，原来的加权曲率条件
\[
\int_0^tW_t(z)(\log f_Z(z))''\,dz<-A_t(0)
\]
恰好等价于 (1.1)。因此中间时间不再是一个独立的曲率积分问题，而是反向风险率导数的比较问题。

还可定义
\[
D_t=
\mathbb E_{\pi_t^-}[\ell_Z]
-\mathbb E_{\pi_t^0}[\ell_Z].
\]
由 (2.1) 及条件 score 恒等式，
\[
\boxed{
D_t=-\frac{h_T'(t)}{h_T(t)}
=-(\log h_T(t))'.
}
\tag{2.4}
\]
故目标也可写成
\[
\boxed{
D_t>\frac{-h_X'(t)}{h_T(t)}.
}
\tag{2.5}
\]

## 3. 一般正次级过程的更新恒等式

设 \(Z\) 是无漂移正次级过程，Laplace 变换为
\[
\widehat f_Z(s)=e^{-\phi_Z(s)},
\qquad
\phi_Z(s)=\int_0^\infty(1-e^{-sx})\Pi_Z(x)\,dx.
\]
定义
\[
g_Z(x)=x\Pi_Z(x).
\]
则
\[
\phi_Z'(s)=
\int_0^\infty e^{-sx}g_Z(x)\,dx
=\mathcal L[g_Z](s).
\]
另一方面，
\[
-\widehat f_Z'(s)
=\mathcal L[x f_Z(x)](s)
=\phi_Z'(s)\widehat f_Z(s).
\]
由 Laplace 卷积定理和唯一性，
\[
\boxed{
x f_Z(x)=(g_Z*f_Z)(x),\qquad x>0.
}
\tag{3.1}
\]
这是正次级过程密度的精确更新方程；其核是 \(g_Z=x\Pi_Z\)，而不是 Lévy 密度本身。

## 4. Robin 延迟谱的显式更新核

沿用
\[
\rho_i=R_i-c,\qquad \sigma_i=L_i-c.
\]
已有精确 Lévy 密度
\[
\Pi_Z(x)=
\sum_{i\ge2}\frac{e^{-\rho_i x}-e^{-\sigma_i x}}{x}.
\]
所以
\[
\boxed{
g_Z(x)=
\sum_{i\ge2}(e^{-\rho_i x}-e^{-\sigma_i x})>0.
}
\tag{4.1}
\]
其 Laplace 变换为
\[
\mathcal L[g_Z](s)
=\sum_{i\ge2}
\left(\frac1{s+\rho_i}-\frac1{s+\sigma_i}\right)
=\phi_Z'(s).
\]

对
\[
\eta_i=\lambda_{i+1}-c
\]
定义两个公共谱更新核
\[
g_X(x)=\sum_{i\ge2}(e^{-\eta_i x}-e^{-\rho_i x}),
\]
\[
g_T(x)=\sum_{i\ge2}(e^{-\eta_i x}-e^{-\sigma_i x}).
\]
在谱排序 \(\eta_i\le\rho_i\le\sigma_i\) 下，三个核均非负，并有
\[
\boxed{g_T=g_X+g_Z.}
\tag{4.2}
\]
相应密度满足
\[
\boxed{
t f_X(t)=(g_X*f_X)(t),
\qquad
t f_T(t)=(g_T*f_T)(t).
}
\tag{4.3}
\]

因此从 \(X_c\) 到 \(T_c=X_c+Z_c\)，更新核的变化恰好是加上一个严格正的 Robin 延迟核 \(g_Z\)。目标 (1.1) 可转化为比较 (4.3) 两个更新方程所产生的反向风险率导数。

## 5. 更新核与半谱计数函数

定义
\[
m(u)=\sum_{i\ge2}\mathbf1_{[\rho_i,\sigma_i)}(u),
\]
以及半谱计数函数
\[
N_R(v)=\#\{i\ge2:R_i\le v\},
\qquad
N_L(v)=\#\{i\ge2:L_i\le v\}.
\]
由于 \(R_i,L_i\) 各自严格递增，且 \(R_i\le L_i\)，固定 \(v=u+c\) 时活跃指标恰为
\[
N_L(v)<i\le N_R(v).
\]
故
\[
\boxed{
m(u)=
\max\{0,N_R(u+c)-N_L(u+c)\}.
}
\tag{5.1}
\]
这说明谱区间重叠数就是两个半谱计数函数之差。

对有限层、正扩散系数的 Sturm–Liouville 半谱，Weyl 计数渐近为
\[
N_R(v)=\frac{\ell_R}{\pi}\sqrt v+O(1),
\qquad
N_L(v)=\frac{\ell_L}{\pi}\sqrt v+O(1).
\]
在 \(\delta_\ell=\ell_R-\ell_L>0\) 的严格光学间隙侧，
\[
\boxed{
m(u)=\frac{\delta_\ell}{\pi}\sqrt u+O(1),
\qquad u\to\infty.
}
\tag{5.2}
\]
这里使用的是 Sturm 振荡理论的谱计数渐近，不涉及逐项相邻谱间距的收敛。

由于
\[
\Pi_Z(x)=\int_0^\infty e^{-ux}m(u)\,du,
\]
由 Watson 型 Laplace 渐近得到
\[
\boxed{
\Pi_Z(x)=
\frac{\delta_\ell}{2\sqrt\pi}x^{-3/2}
+O(x^{-1}),
\qquad x\downarrow0.
}
\tag{5.3}
\]
于是
\[
\boxed{
g_Z(x)=
\frac{\delta_\ell}{2\sqrt{\pi x}}+O(1),
\qquad x\downarrow0.
}
\tag{5.4}
\]
这从谱计数直接重获了小时间光学长度控制的高频主项。

在当前静态打平基准 \(\rho_2<\sigma_2<\rho_3\) 下，
\[
g_Z(x)\sim e^{-\rho_2x},\qquad x\to\infty.
\tag{5.5}
\]
更一般地，晚时指数由第一个严格非零谱区间的左端点决定。

## 6. 与 Bondesson、GGC 及广义单峰理论的关系

由于 \(\Pi_Z\) 严格完全单调，\(Z_c\) 属于具有完全单调 Lévy 密度的 Bondesson 类。可参见：

- L. Bondesson, *Classes of Infinitely Divisible Distributions and Densities*, Z. Wahrscheinlichkeitstheorie verw. Gebiete 57 (1981), 39–71, [DOI](https://doi.org/10.1007/BF00533713).

GGC 是其中更窄的 Thorin 类，其 Bernstein 指数导数 \(\phi_Z'\) 必须是 Stieltjes 函数。这里
\[
\phi_Z'(s)=\int_0^\infty\frac{m(u)}{(s+u)^2}\,du.
\]
分部积分给出形式上的 Stieltjes 表示
\[
\phi_Z'(s)=\int_0^\infty\frac{dm(u)}{s+u}.
\]
因此，若 \(m\) 非递减，便得到正 Thorin 测度。反过来，Stieltjes 表示的唯一性要求相应的 \(dm\) 为非负测度。

但当前静态打平基准已有严格区间证书
\[
m(u_+)\ge13,\qquad m(u_-)\le11,
\qquad u_+<u_-.
\]
所以 \(m\) 不是非递减函数，\(dm\) 含有负向跳跃，\(\phi_Z'\) 不是 Stieltjes 函数。由此
\[
\boxed{Z_c\text{ 在该基准下不是 GGC。}}
\]
Bondesson 文献中 GGC 密度的单峰性定理因此不能直接用于本模型。

Bridaa–Jedidi–Sendov (2024), *Generalized Unimodality and Subordinators*, [DOI](https://doi.org/10.1007/s10959-023-01242-z)，研究的是 Bernstein 指数的 \(BF_s\) 类及相应的 star-unimodality。其 Lévy 密度要求具有 \(p_s(x)x^{-s-1}\) 形式，其中 \(p_s\) 非递减。当前 \(\Pi_Z(x)\) 具有指数衰减尾部，\(x^{s+1}\Pi_Z(x)\to0\)，不可能是一个处处非递减且严格为正的 \(p_s\)。因此该类也不直接覆盖当前 \(Z_c\)。

## 7. 小时间、晚时间与中间时间

由已有 Robin 谱渐近，令
\[
d_X=\ell_L,\qquad d_T=\ell_L+\delta_\ell=\ell_R.
\]
则
\[
h_X'(t)\sim-\frac{d_X^2}{2t^3},
\qquad
h_T'(t)\sim-\frac{d_T^2}{2t^3}.
\]
所以
\[
\boxed{
h_X'(t)-h_T'(t)
\sim
\frac{d_T^2-d_X^2}{2t^3}>0,
\qquad t\downarrow0.
}
\tag{7.1}
\]

晚时间沿用最慢极点留数：
\[
h_X'(t)\sim-\eta B_Xe^{-\eta t},
\qquad
h_T'(t)\sim-\eta B_Te^{-\eta t},
\]
其中
\[
\frac{B_T}{B_X}
=\mathbb E[e^{\eta Z_c}]
=:M>1.
\]
因此
\[
\boxed{
h_X'(t)-h_T'(t)
\sim
\eta B_X(M-1)e^{-\eta t}>0.
}
\tag{7.2}
\]
归一化裕量为
\[
\boxed{
\frac{h_X'(t)-h_T'(t)}{h_T(t)}
\longrightarrow
\eta\left(1-\frac1M\right)>0.
}
\tag{7.3}
\]

故固定非退化 Robin 几何的小时间与充分晚时间均已关闭。仍然需要证明的是
\[
\boxed{
h_X'(t)-h_T'(t)>0
}
\]
在整个有限中间时间区间成立，尤其是 \(c-a\downarrow0\) 时的参数统一版本。

## 8. 下一步的精确任务

现在有两个互补的描述：

1. 条件概率描述：
   \[
   h_X'-h_T'
   =
   \mathbb E_{\pi_t}A_t-A_t(0)
   -\operatorname{Var}_{\pi_t}(q_t).
   \]
2. 谱更新描述：
   \[
   t f_X=g_X*f_X,\qquad
   t f_T=(g_X+g_Z)*f_T.
   \]

下一步应研究一个模型专属的比较命题：

**当更新核由 \(g_X\) 增加为 \(g_X+g_Z\)，且 \(g_Z\) 具有 Robin 谱计数结构时，反向风险率导数是否严格下降？**

一般正次级过程并不自动提供这一结论，所以必须利用 (5.1)–(5.4) 的具体谱结构。目标仍是证明
\[
h_T'(t)<h_X'(t),
\]
然后再接入完整 \((\log r)''\)。

本轮没有重新计算谱根或进行数值扫描；所有小时间、晚时间留数和既有谱证书均直接沿用。
