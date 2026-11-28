# Robin 延迟谱乘积的 Lévy 表示与 bell-shaped 密度判据

## 1. 对象

沿用指数倾斜后的延迟变量
\[
\widehat\nu_c(s)=
\prod_{i\ge2}\frac{1+s/\sigma_i}{1+s/\rho_i},
\qquad
\rho_i=R_i-c,\quad \sigma_i=L_i-c,
\]
其中
\[
0<\rho_i\le \sigma_i,\qquad
\rho_i,\sigma_i\asymp i^2.
\]
严格光学长度差给出
\[
\prod_{i\ge2}\frac{\rho_i}{\sigma_i}=0,
\]
故 \(Z_c\) 无零点原子。这里不重新研究其晚时留数或小时间鞍点。

## 2. 精确 Lévy–Khintchine 表示

逐项使用
\[
\log\frac{1+s/\rho_i}{1+s/\sigma_i}
=
\int_{\rho_i}^{\sigma_i}\frac{s}{u(u+s)}\,du
=
\int_0^\infty(1-e^{-sx})
\left(\int_{\rho_i}^{\sigma_i}e^{-ux}\,du\right)dx,
\]
得到
\[
\boxed{
-\log\widehat\nu_c(s)
=
\int_0^\infty(1-e^{-sx})\Pi_c(x)\,dx,
}
\]
其中
\[
\boxed{
\Pi_c(x)=
\sum_{i\ge2}\int_{\rho_i}^{\sigma_i}e^{-ux}\,du
=
\sum_{i\ge2}\frac{e^{-\rho_i x}-e^{-\sigma_i x}}{x}.
}
\tag{2.1}
\]
因此 \(Z_c\) 是无漂移次级过程在单位时间的分布。

对每个 \(n\ge0\)，
\[
(-1)^n\Pi_c^{(n)}(x)
=
\sum_{i\ge2}\int_{\rho_i}^{\sigma_i}u^ne^{-ux}\,du>0.
\]
所以
\[
\boxed{\Pi_c\text{ 是严格完全单调的 Lévy 密度。}}
\]
这将 \(Z_c\) 归入具有完全单调 Lévy 密度的 Bondesson 类。

定义谱区间重叠计数
\[
\boxed{
m(u)=\sum_{i\ge2}\mathbf1_{[\rho_i,\sigma_i)}(u).
}
\]
则
\[
\Pi_c(x)=\int_0^\infty e^{-ux}m(u)\,du.
\]
由 \(\rho_i\asymp i^2\)，
\[
\int_1^\infty\frac{m(u)}{u^3}\,du
=
\sum_{i\ge2}\int_{\max(1,\rho_i)}^{\sigma_i}u^{-3}\,du
\le \frac12\sum_{i\ge2}\rho_i^{-2}<\infty.
\tag{2.2}
\]
同时
\[
\int_1^\infty\frac{m(u)}{u^2}\,du
\le\sum_{i\ge2}\rho_i^{-1}<\infty.
\tag{2.3}
\]

## 3. Kwaśnicki–Simon bell-shape 判据

取特征函数
\[
\Phi(\xi)=\mathbb E e^{-i\xi Z_c}=\widehat\nu_c(i\xi).
\]
由 Lévy 表示，
\[
\log\Phi(\xi)
=
\int_0^\infty
\left(\frac1{u+i\xi}-\frac1u\right)m(u)\,du.
\]
将 \(u>1\) 的一阶项分离：
\[
\begin{aligned}
\log\Phi(\xi)
={}&-ib\xi\\
&+\int_0^\infty
\left[
\frac1{u+i\xi}
-\left(\frac1u-\frac{i\xi}{u^2}\right)\mathbf1_{\{u>1\}}
\right]m(u)\,du,
\end{aligned}
\]
其中
\[
b=\int_1^\infty\frac{m(u)}{u^2}\,du<\infty.
\]
这正是 Kwaśnicki 的 Stieltjes 型 Fourier 表示，且其 level-crossing 函数为
\[
\varphi(u)=
\begin{cases}
m(u),&u>0,\\
0,&u<0.
\end{cases}
\]
正则条件自动成立：\(|\Phi(\xi)|\le1\)，故
\[
\int_{-1}^1\operatorname{Re}\Phi(\xi)\,d\xi<\infty,
\qquad
\lim_{\xi\to0}\xi\operatorname{Im}\Phi(\xi)=0.
\]
可积条件由 (2.2) 给出。

Kwaśnicki 的定理指出：若对每个整数 \(k\)，\(\varphi-k\) 至多改变一次符号，且 \(k=0\) 时满足正负半轴的规定符号，则对应函数 weakly bell-shaped；若该密度光滑，则它是 bell-shaped。参见：

- M. Kwaśnicki, *A new class of bell-shaped functions*, Transactions of the AMS 373 (2020), 349–382, [DOI](https://doi.org/10.1090/tran/7825).
- M. Kwaśnicki and T. Simon, *Characterisation of the class of bell-shaped functions*, Mathematische Zeitschrift 301 (2022), 3207–3230, [DOI](https://doi.org/10.1007/s00209-022-02997-7).

在本模型中，\(m(u)\ge0\) 且 \(m(u)=0\) 对 \(u<\rho_2\)，所以 \(k=0\) 条件自动满足。真正需要检查的是整数水平 \(k\ge1\)。

## 4. 整数水平穿越的精确谱条件

因为 \((\rho_i)\) 与 \((\sigma_i)\) 均严格递增，所有在给定 \(u\) 处活跃的谱区间指标构成连续整数段。于是
\[
\boxed{
m(u)\ge k+1
\iff
\exists i\ge2:\quad \rho_{i+k}\le u<\sigma_i.
}
\tag{4.1}
\]
因此 \(m(u)-k\) 的 level-crossing 条件等价于以下谱区间条件：

\[
\boxed{
\mathcal I_k:=\{i\ge2:\rho_{i+k}<\sigma_i\}
\text{ 对每个 }k\ge1\text{ 都是一个上尾区间（或空集）。}
}
\tag{4.2}
\]

这里“上尾区间”允许端点处相等；相等只对应两个谱端点同时到达，不产生负的水平穿越。

解释：\(\rho_{i+k}<\sigma_i\) 表示第 \(i\) 至第 \(i+k\) 个谱区间曾同时覆盖某个频率。如果该关系在某个指标之后失效，则在后续 \(\sigma_i\) 处会有少于 \(k\) 个区间仍然活跃，导致 \(m-k\) 从正值回到负值。反之，若 (4.2) 成立，一旦 \(m\) 超过 \(k\)，后续每次区间退出都至少保留 \(k\) 个活跃区间，故不会再次跌到 \(k\) 以下。

## 5. 为什么该条件可化为有限低模检查

令
\[
\Delta_i^L=\sqrt{L_{i+1}}-\sqrt{L_i},\qquad
\Delta_j^R=\sqrt{R_{j+1}}-\sqrt{R_j}.
\]
对有限层、正扩散系数的 Sturm–Liouville 半谱，
\[
\Delta_i^L\longrightarrow\frac{\pi}{\ell_L},
\qquad
\Delta_j^R\longrightarrow\frac{\pi}{\ell_R}.
\]
在本研究的严格光学间隙侧，
\[
\ell_L<\ell_R,
\]
故
\[
\frac{\pi}{\ell_L}>\frac{\pi}{\ell_R}.
\]
于是存在一个有限整数 \(I_0\)，使得
\[
\inf_{i\ge I_0}\Delta_i^L
>
\sup_{j\ge I_0}\Delta_j^R.
\tag{5.1}
\]

对固定 \(k\ge1\)，定义
\[
d_{i,k}=\sqrt{L_i}-\sqrt{R_{i+k}}.
\]
当 \(i\ge I_0\) 时，
\[
d_{i+1,k}-d_{i,k}
=\Delta_i^L-\Delta_{i+k}^R>0.
\tag{5.2}
\]
所以对每个 \(k\)，\(d_{i,k}\) 在 \(i\ge I_0\) 上严格递增，并最终趋于 \(+\infty\)。又因 \(R_{i+k}\to\infty\)，对所有 \(i<I_0\)，当 \(k\) 足够大时 \(d_{i,k}<0\)。

定义有限阈值
\[
K_0=
\max_{2\le i<I_0}
\max\{k\ge1:\sqrt{R_{i+k}}<\sqrt{L_i}\},
\]
若集合为空则令 \(K_0=0\)。那么：

- 对 \(k>K_0\)，所有 \(i<I_0\) 均无重叠；从 \(I_0\) 起 \(d_{i,k}\) 严格递增，故 \(\mathcal I_k\) 自动为上尾区间。
- 对 \(1\le k\le K_0\)，只需检查有限个低模指标 \(2\le i\le I_0\)，确认 \(\mathcal I_k\) 在该有限前缀中没有“真—假”回落；之后由 (5.2) 自动保持。

因此，**Kwaśnicki–Simon 定理的全部整数水平条件，可由高频谱间距定理加有限个低模区间检查完成。**

这是一项重要的简化：不需要对无穷多个 \(k\) 和无穷多个谱模逐一进行独立认证。

## 6. 对 \(f_Z\) 的结论及其对当前证明的价值

若有限低模检查通过，并采用已有正光学距离的 sectorial 高频估计
\[
|\widehat\nu_c(i\xi)|
\le C e^{-(\ell_R-\ell_L)\sqrt{|\xi|/2}}
\]
（大 \(|\xi|\)），则 \(\Phi\) 的任意多项式加权均可积，Fourier 反演给出 \(f_Z\in C^\infty(\mathbb R)\)。因此 Kwaśnicki–Simon 判据推出
\[
\boxed{f_Z\text{ 是严格 bell-shaped 密度。}}
\]

特别地：

1. \(f_Z'\) 恰有一个符号变化，故 \(f_Z\) 严格单峰。
2. \(f_Z''\) 恰有两个符号变化，故密度有且仅有两个拐点。
3. 在 \(f_Z''<0\) 的中央区间，
   \[
   (\log f_Z)''
   =\frac{f_Z''}{f_Z}
    -\left(\frac{f_Z'}{f_Z}\right)^2<0.
   \]
4. 已有谱留数结果表明，在特定 Robin 谱交错下，\((\log f_Z)''\) 的晚时尾部严格为正。

因此，bell-shape 结构将延迟密度的普通曲率符号限制在两个外侧凸区间和一个中央凹区间中。它没有直接证明 \((\log f_Z)''\) 只有一次变号，也没有直接推出 \(r(t)\) 单峰；这两点不能从 bell-shaped 定义中跳过。

## 7. 与最终曲率问题的连接

仍然使用已经证明的恒等式
\[
h_T'-h_X'
=
A_t(0)+
\int_0^tW_t(z)(\log f_Z(z))''\,dz,
\]
其中
\[
W_t'(z)=p_t(z)[h_T-h_X(t-z)]
\]
且 \(W_t\) 严格单峰。

本轮新结果提供了一个可用的分区：

- \(f_Z''<0\) 的中央区间自动贡献负的对数曲率积分；
- 晚时 \(f_Z''>0\) 区间需要继续使用已有谱留数上界；
- 小时间外侧区间由高频鞍点估计处理；
- 只需进一步定位 \(W_t\) 的峰值 \(y_t\) 与两个密度拐点的相对位置，便可尝试将积分分成有明确符号的区域。

这是一条比全局对数凹性更弱、也更符合已知 Robin 谱结构的路线。

## 8. 本轮结论与待完成事项

**本轮已解析建立：**

1. \(Z_c\) 的精确 Lévy–Khintchine 表示。
2. Lévy 密度的严格完全单调性。
3. Kwaśnicki–Simon bell-shape 定理在本谱乘积中的精确 level-crossing 函数。
4. 整数水平条件与 Robin 谱区间重叠条件 (4.2) 的等价关系。
5. 高频谱间距将无限水平条件压缩为有限低模检查的充分条件。

**尚未完成：**

- 对当前整个 Robin 参数域证明有限低模重叠检查通过；
- 利用 bell-shape 得到 \((\log f_Z)''\) 的单次变号或足够强的尾部积分界；
- 证明完整的 \((\log r(t))''<0\)。

本轮没有重新计算谱根，也没有把条件性 bell-shape 结论误写成整个参数域已通过的定理。
