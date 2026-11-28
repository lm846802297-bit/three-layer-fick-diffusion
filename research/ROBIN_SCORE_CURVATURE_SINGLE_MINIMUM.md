# Robin 延迟密度 score 的唯一极小点：谱留数表示与变差障碍

## 1. 本阶段的精确命题

令
\[
\Phi_Z(s)=\prod_{j\ge2}\frac{1+s/\sigma_j}{1+s/\rho_j},
\qquad
0<\rho_j<\sigma_j,\qquad
\rho_j,\sigma_j\asymp j^2,
\]
其中 \(\rho_j=R_j-c\)、\(\sigma_j=L_j-c\)，并假设谱根简单且没有零极点相消。记其正密度为 \(f_Z\)，score 与 score 曲率为
\[
\ell_Z(t)=(\log f_Z(t))',
\qquad
G_Z(t)=\ell_Z'(t)
=(\log f_Z(t))''.
\]

当前严格已知的两端符号为
\[
G_Z(t)\sim-\frac{\delta_\ell^2}{2t^3}<0,
\qquad t\downarrow0,
\tag{1.1}
\]
其中 \(\delta_\ell=\ell_R-\ell_L>0\)，以及在
\[
\rho_2<\sigma_2<\rho_3,\qquad \sigma_3>\rho_3
\]
等使前两个留数为正的谱条件下，
\[
G_Z(t)\sim
\frac{A_3}{A_2}(\rho_3-\rho_2)^2
e^{-(\rho_3-\rho_2)t}>0,
\qquad t\to\infty.
\tag{1.2}
\]
因此至少存在一个 \(z_*>0\) 使 \(G_Z(z_*)=0\)。

本阶段要证明的是：这个零点唯一，且 \(G_Z\) 从负变正。全局对数凹性和全局对数凸性均不是目标。

## 2. 曲率分子有精确的双谱留数展开

在 Laplace 反演的留数级数及其前两阶导数于每个紧子区间 \(t\ge\varepsilon>0\) 正常收敛的条件下，
\[
f_Z(t)=\sum_{j\ge2}A_j e^{-\rho_jt},
\]
其中
\[
A_j=
\rho_j\left(1-\frac{\rho_j}{\sigma_j}\right)
\prod_{\substack{i\ge2\\i\ne j}}
\frac{1-\rho_j/\sigma_i}{1-\rho_j/\rho_i}.
\tag{2.1}
\]
定义
\[
Q_Z(t)=f_Z(t)f_Z''(t)-f_Z'(t)^2.
\]
逐项相乘并将 \((i,j)\)、\((j,i)\) 配对，得到
\[
\boxed{
Q_Z(t)=
\sum_{2\le i<j}
A_iA_j(\rho_j-\rho_i)^2
e^{-(\rho_i+\rho_j)t}.
}
\tag{2.2}
\]
由于 \(f_Z(t)>0\)，
\[
\boxed{
\operatorname{sgn}G_Z(t)=\operatorname{sgn}Q_Z(t).
}
\tag{2.3}
\]
所以唯一 score 极小点问题等价于证明双谱指数级数 \(Q_Z\) 在正半轴恰有一个零点。

式 (2.2) 是本阶段的核心表示：它把曲率问题从对密度求导，转化成由两两谱率之和 \(\rho_i+\rho_j\) 排序的有符号指数级数。

## 3. 留数符号由谱区间重叠数精确决定

定义
\[
m(u)=\sum_{i\ge2}\mathbf1_{[\rho_i,\sigma_i)}(u).
\]
在不存在端点重合的情形，令
\[
m_j=m(\rho_j+)
=1+\#\{i<j:\sigma_i>\rho_j\}.
\]
在 (2.1) 的留数乘积中，\(i>j\) 的因子均为正；对 \(i<j\)，分母均为负，而分子为正恰好当 \(\sigma_i>\rho_j\)。因此
\[
\boxed{
\operatorname{sgn}A_j=(-1)^{m_j-1}.
}
\tag{3.1}
\]
这给出了一个直接联系：谱区间在 \(\rho_j\) 右侧的重叠数奇偶性，决定第 \(j\) 个极点留数的符号。

### 3.1 严格光学长度差迫使留数符号无限次改变

令
\[
k_j=\#\{i\ge2:\rho_j<\sigma_i<\rho_{j+1}\}.
\]
在端点不重合时，从 \(\rho_j+\) 移动到 \(\rho_{j+1}+\)，一个新谱区间进入，同时有 \(k_j\) 个旧区间退出，所以
\[
\boxed{m_{j+1}=m_j+1-k_j.}
\tag{3.2}
\]

由两侧 Weyl 计数，
\[
m_j
=
\#\{i\ge2:R_i\le R_j\}
-\#\{i\ge2:L_i\le R_j\}
=
\frac{\ell_R-\ell_L}{\ell_R}\,j+O(1).
\tag{3.3}
\]
故 \(m_j\to\infty\)。

如果从某个指标起每个区间 \((\rho_j,\rho_{j+1})\) 都至少包含一个 \(\sigma_i\)，即 \(k_j\ge1\)，那么 (3.2) 会使 \(m_j\) 从该指标起单调不增，与 (3.3) 矛盾。因此
\[
\boxed{k_j=0\quad\text{对无穷多个 }j.}
\tag{3.4}
\]
在这些间隙中 \(m_{j+1}=m_j+1\)，于是由 (3.1)
\[
\operatorname{sgn}A_{j+1}=-\operatorname{sgn}A_j.
\]
所以
\[
\boxed{\text{Robin 延迟谱的留数 }A_j
\text{ 在高阶谱尾中必然发生无穷多次符号改变。}}
\tag{3.5}
\]

这项结论只依赖严格光学长度差和两侧 Weyl 计数，不需要重新计算谱根。

## 4. 一个可用但目前不能直接满足的单零点判据

将 (2.2) 写成有符号离散测度的 Laplace 变换：
\[
Q_Z(t)=\int_0^\infty e^{-\lambda t}\,\mu_Q(d\lambda),
\]
其中
\[
\mu_Q=
\sum_{i<j}
A_iA_j(\rho_j-\rho_i)^2
\delta_{\rho_i+\rho_j}.
\tag{4.1}
\]

**单变差充分条件。** 若存在 \(\lambda_*>0\)，使 \(\mu_Q\) 的正部分支撑在 \((0,\lambda_*]\)，负部分支撑在 \([\lambda_*,\infty)\)，并且两部分均非零，则 \(Q_Z(t)\) 至多有一个零点。若再结合 (1.1)–(1.2)，该零点恰好一个，且符号由负变正。

证明：写
\[
P(t)=\int e^{-\lambda t}\mu_Q^+(d\lambda),
\qquad
N(t)=\int e^{-\lambda t}\mu_Q^-(d\lambda).
\]
支撑次序给出
\[
\frac{d}{dt}\log\frac{P(t)}{N(t)}
=
\mathbb E_{N,t}[\lambda]
-\mathbb E_{P,t}[\lambda]\ge0.
\]
因此 \(P/N\) 单调，方程 \(P=N\) 至多有一个解。结合两端异号即得唯一零点。

**当前障碍。** (3.5) 表明 \(A_j\) 不可能在某个高阶指标之后保持固定符号。因此，不能通过“所有高阶留数同号”来证明 (4.1) 的单变差条件。固定 \(i=2\) 的子族
\[
A_2A_j(\rho_j-\rho_2)^2
\delta_{\rho_2+\rho_j}
\]
本身就继承了无限次留数符号变化。其它谱对可能与它们在指数率上交错或发生聚合，故目前还不能据此断言完整 \(\mu_Q\) 必然有无限次变差；但简单的单变差证明路线已经没有依据。

## 5. 与已有概率理论的关系

1. Sumita–Masuda (1987) 研究 Laplace 变换仅有负实零点与极点的密度类 \(\Omega_+\)，其中包括若干独立指数变量的和与混合；其摘要还涉及 modified renewal densities 的完全单调性和单峰性。当前 Robin 乘积属于相邻的负实零点/极点谱框架，但该文公开摘要没有给出 \((\log f_Z)''\) 的唯一变号定理。
2. Hansen (1988), Theorems 3–4，给出正无穷可分密度全局对数凹/对数凸的 Lévy 密度条件。这里 \(\Pi_Z\) 严格完全单调，因而严格对数凸；并且
   \[
   \Pi_Z(t)\sim
   \frac{\delta_\ell}{2\sqrt\pi}t^{-3/2}\to\infty
   \quad(t\downarrow0).
   \]
   所以 Hansen 的全局对数凹/凸判据均不覆盖当前密度。它们没有解决局部曲率只变号一次的问题。
3. Keilson (1971) 的 passage-time 对数凹/凸结果针对特定扩散与生灭过程的 passage-time 密度；不能直接推出当前 Robin 延迟谱密度的 score 曲率唯一变号。

参考文献：

- U. Sumita and Y. Masuda, *Classes of probability density functions having Laplace transforms with negative zeros and poles*, Advances in Applied Probability 19 (1987), 632–651. [DOI](https://doi.org/10.2307/1427410).
- B. G. Hansen, *On Log-Concave and Log-Convex Infinitely Divisible Sequences and Densities*, Annals of Probability 16 (1988), 1832–1839. [DOI](https://doi.org/10.1214/aop/1176991600).
- J. Keilson, *Log-concavity and log-convexity in passage time densities of diffusion and birth-death processes*, Journal of Applied Probability 8 (1971), 391–398. [DOI](https://doi.org/10.2307/3211909).

## 6. 当前结论与下一步

本轮已经严格得到：

1. score 曲率的双谱留数表达式 (2.2)；
2. 单个留数的符号由谱重叠计数奇偶性决定，见 (3.1)；
3. 严格光学长度差迫使高阶留数发生无穷多次符号改变，见 (3.5)；
4. 单变差 Laplace 判据足以推出唯一曲率转折点，但当前谱结构不能通过固定高阶留数符号来满足它。

因此，唯一 score 极小点仍未证明，也没有被否定。真正需要的新估计是控制 (2.2) 中高阶成对谱项的**聚合贡献**，而不是逐个控制留数符号。

下一步优先研究截断分解
\[
Q_Z(t)=Q_{\le N}(t)+Q_{>N}(t),
\]
其中低模项由有限个留数精确控制，高模项由 Robin 高频传递矩阵与 Weyl 余项给出统一上界。目标是在一个固定中间时间窗口内证明高模聚合项不产生额外零点。全局对数凹性、bell-shaped 性质和重新计算已有谱根均不在本阶段任务内。

本轮只使用已有谱渐近与留数结论，没有启动新的谱根计算或数值扫描。
