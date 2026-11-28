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


## 9. Score-renewal 方程与曲率权重的历史输运

本节继续使用
\[
x f_Z(x)=(g_Z*f_Z)(x),
\qquad
g_Z(x)=x\Pi_Z(x).
\]
对固定 \(x>0\)，定义后验随机变量 \(Y_x\)：
\[
\boxed{
\Pr(Y_x\in dy)
=
\frac{g_Z(y)f_Z(x-y)}
{x f_Z(x)}\,\mathbf1_{(0,x)}(y)\,dy.
}
\tag{9.1}
\]
更新方程保证它是概率密度。令 \(V_x=x-Y_x\)，并记
\[
\ell_Z(x)=(\log f_Z(x))'.
\]

由于 Robin 正光学间隙给出 \(f_Z(0+)=0\)，且 \(g_Z\) 在原点局部可积，对卷积恒等式求导可得
\[
\frac{(g_Z*f_Z)'(x)}{(g_Z*f_Z)(x)}
=
\mathbb E[\ell_Z(V_x)].
\]
另一方面 \(g_Z*f_Z=x f_Z\)，故
\[
\boxed{
\ell_Z(x)+\frac1x
=
\mathbb E[\ell_Z(V_x)].
}
\tag{9.2}
\]

对 (9.2) 关于 \(x\) 求导。后验密度 (9.1) 的对数导数为
\[
\partial_x\log p_x(y)
=
\ell_Z(x-y)-\ell_Z(x)-\frac1x,
\]
其均值为零。于是
\[
\frac{d}{dx}\mathbb E[\ell_Z(V_x)]
=
\mathbb E[\ell_Z'(V_x)]
+
\operatorname{Var}(\ell_Z(V_x)).
\]
得到精确 score-renewal 恒等式
\[
\boxed{
(\log f_Z)''(x)
=
\mathbb E[(\log f_Z)''(V_x)]
+
\operatorname{Var}(\ell_Z(V_x))
+\frac1{x^2}.
}
\tag{9.3}
\]
在当前 Robin 小时间鞍点渐近下，后验 score 的一阶、二阶矩在每个固定 \(x>0\) 有限，因此上述微分与交换积分合法。

### 9.1 将曲率权重向早时间输运

把 (9.3) 乘以 \(W_t(x)\)，在 \(0<x<t\) 积分。定义
\[
\boxed{
\widetilde W_t(v)
=
f_Z(v)\int_v^t
\frac{W_t(x)g_Z(x-v)}
{x f_Z(x)}\,dx.
}
\tag{9.4}
\]
由 Fubini 定理，
\[
\boxed{
\begin{aligned}
\int_0^t W_t(x)(\log f_Z)''(x)\,dx
={}&
\int_0^t\widetilde W_t(v)(\log f_Z)''(v)\,dv\\
&+
\int_0^tW_t(x)
\left[
\operatorname{Var}(\ell_Z(V_x))+\frac1{x^2}
\right]dx.
\end{aligned}
}
\tag{9.5}
\]
并且
\[
\boxed{
\int_0^t\widetilde W_t(v)\,dv
=
\int_0^tW_t(x)\,dx.
}
\tag{9.6}
\]
其概率解释是：先按 \(W_t(x)\,dx\) 选择当前时间 \(x\)，再按更新后验选择 \(V_x=x-Y_x\)。由于 \(0<V_x<x\)，归一化后的 \(\widetilde W_t\) 相对于 \(W_t\) 严格向早时间输运；特别地，对任意递增函数 \(\psi\)，
\[
\frac{\int_0^t\widetilde W_t(v)\psi(v)\,dv}
{\int_0^tW_t(x)\,dx}
\le
\frac{\int_0^tW_t(x)\psi(x)\,dx}
{\int_0^tW_t(x)\,dx}.
\tag{9.7}
\]

式 (9.5) 把原曲率积分拆成历史曲率的输运贡献与两个严格正的更新补偿项。它没有单独给出目标不等式的符号，因为 \((\log f_Z)''\) 并非全局递增；其用途是将晚时正曲率如何由早时负曲率及更新方差补偿，转化为一个具有显式后验核的比较问题。

## 10. 当前证明状态

已得到的精确对象为：

1. 原加权曲率积分等于 \(h_T'(t)\)；
2. \(Z_c\) 的 Lévy 更新核 \(g_Z=x\Pi_Z\)；
3. 两个累计核满足 \(g_T=g_X+g_Z\)；
4. 更新核由半谱计数差 \(m(u)=N_R-N_L\) 决定；
5. score 满足历史更新恒等式 (9.3)；
6. 单峰权重 \(W_t\) 在更新后验下输运为同质量、偏向早时间的 \(\widetilde W_t\)。

仍需证明的是更新后验所产生的正方差补偿，如何与 \(\widetilde W_t\) 对早时负曲率的加权作用共同给出
\[
h_X'(t)-h_T'(t)>0
\]
的中间时间统一下界。本节不将该比较写成已证明结论。


## 11. 中间时间的单峰权重分割与可检验尾项界

令
\[
\ell_Z=(\log f_Z)',\qquad
\mathcal G(t)=h_X'(t)-h_T'(t).
\]
精确公式为
\[
\mathcal G(t)=A_t(0)+\int_0^tW_t(z)\ell_Z'(z)\,dz.
\tag{11.1}
\]
若能证明 Robin 延迟密度的 score \(\ell_Z\) 有唯一极小点 \(z_*\)，即
\[
\ell_Z'(z)<0\quad (z<z_*),\qquad
\ell_Z'(z)>0\quad (z>z_*),
\tag{11.2}
\]
则对 \(t>z_*\)，定义
\[
L_t=\int_0^{z_*}W_t(z)(-\ell_Z'(z))\,dz,
\qquad
R_t=\int_{z_*}^tW_t(z)\ell_Z'(z)\,dz.
\]
于是
\[
\boxed{\mathcal G(t)=A_t(0)+L_t-R_t.}
\tag{11.3}
\]
所以中间时间的充分且必要条件是
\[
\boxed{L_t>R_t-A_t(0).}
\tag{11.4}
\]
(11.2) 是一个独立的辅助命题，不能由晚时对数凸性单独推出；在本模型中仍须由 Robin 谱结构证明，或者用不依赖唯一曲率转折点的直接估计替代。

### 11.1 单峰权重给出的严格尾项上界

已有
\[
W_t'(z)=p_t^-(z)[h_T(t)-h_X(t-z)]
\]
且 \(h_X(t-z)\) 关于 \(z\) 严格递增。因此 \(W_t\) 在唯一的 \(y_t\) 处达到最大值：
\[
W_t(z)\le W_t(y_t),\qquad 0<z<t.
\]
若 (11.2) 成立，则
\[
\boxed{
R_t
\le
W_t(y_t)\int_{z_*}^t\ell_Z'(z)\,dz
=
W_t(y_t)[\ell_Z(t)-\ell_Z(z_*)].
}
\tag{11.5}
\]
这给出一个不需要估计整条晚时曲率曲线的尾项上界。由条件 CDF 差表示，
\[
W_t(y_t)
=h_T(t)\left[
F_{\pi_t^-}(y_t)-F_{\pi_t^0}(y_t)
\right].
\tag{11.6}
\]
故可将 (11.5) 写成完全概率化的上界：
\[
R_t\le h_T(t)\,
\bigl[F_{\pi_t^-}(y_t)-F_{\pi_t^0}(y_t)\bigr]\,
[\ell_Z(t)-\ell_Z(z_*)].
\tag{11.7}
\]

因此，一个明确的充分条件是
\[
\boxed{
L_t>
h_T(t)\,
\bigl[F_{\pi_t^-}(y_t)-F_{\pi_t^0}(y_t)\bigr]\,
[\ell_Z(t)-\ell_Z(z_*)]
-A_t(0).
}
\tag{11.8}
\]
它把正曲率尾项压缩为三个可独立估计的因子：反向风险率尺度、条件 CDF 最大差、score 回升幅度。

### 11.2 早时负曲率的可认证下界

对任意 \(0<a<b\le\min\{z_*,y_t\}\)，因 \(W_t\) 在 \((0,y_t)\) 严格递增，
\[
L_t\ge
W_t(a)\int_a^b[-\ell_Z'(z)]\,dz
=
W_t(a)[\ell_Z(a)-\ell_Z(b)].
\tag{11.9}
\]
结合 (11.8)，得到一个便于分段认证的充分判据：
\[
\boxed{
\begin{aligned}
&W_t(a)[\ell_Z(a)-\ell_Z(b)]\\
&\quad>
W_t(y_t)[\ell_Z(t)-\ell_Z(z_*)]-A_t(0).
\end{aligned}
}
\tag{11.10}
\]
其中 \(W_t(a)\) 和 \(W_t(y_t)\) 均可由两个条件分布的 CDF 差直接计算或估计。

这不是对所有 \(t\) 已完成的证明；它给出了可用于解析区间分割的明确验收条件。要关闭中间时间，需要证明 (11.10) 在覆盖区间成立，或改进其上下界。

### 11.3 文献边界与本模型所需的新内容

已有可靠性文献中，DRHR 与 log-concavity of the CDF 是标准等价关系；Sangüesa, Badía and Cha (2014) 的 Proposition 3.1(b) 研究独立正增量磨损过程下 DRHR 的保持，其假设包括过程处于 IPDI 类以及相应的 hazard-rate stochastic ordering。该定理不是对任意正次级过程卷积密度 score 的单谷定理，也没有给出 (11.8) 的定量常数。

Alimohammadi, Balakrishnan and Simon (2024) 的 Theorem 2.1 给出 log-concave 函数的积分刻画，并用于 failure-rate 单调性；它没有直接控制当前 \(W_t\) 与 \((\log f_Z)''\) 的加权竞争。

所以现成文献可以提供随机序与积分不等式工具，但当前 Robin 延迟谱需要额外证明的仍是 (11.2) 或直接的 (11.10)，以及其中的参数统一常数。

## 12. 本阶段验收条件

当前阶段只按以下标准验收，不将后续 \(r(t)\) 结论提前并入：

1. 对固定 Robin 几何，证明或替代 (11.2) 的 score 曲率结构；
2. 在整个剩余中间时间区间，解析证明 (11.4)，优先使用 (11.10)；
3. 对首模逆序趋于打平的参数族，给出所需常数的统一控制；
4. 得到严格结论 \(h_T'(t)<h_X'(t)\)。

达到第 4 项后，再将该结论代回完整 \((\log r(t))''\) 恒等式，处理其中其余谱项并证明 \(r(t)\) 单峰。当前不把 \(h_T'<h_X'\) 单独等同于完整的 \(r\)-单峰定理。
