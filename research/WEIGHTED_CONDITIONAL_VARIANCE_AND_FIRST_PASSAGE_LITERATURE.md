# 条件方差的加权控制、首次通过时间文献边界与近打平晚时渐近

## 1. 本轮目标

沿用首模唯一逆序条件
\[
a=L_1<c=\lambda_2<b=R_1,
\qquad
L_i\ge R_i\ge\lambda_{i+1}\quad(i\ge2).
\]
令 \(X_c,Z_c\) 独立，\(T_c=X_c+Z_c\)，其 CDF 分别为 \(F_X,F_T\)，反向风险率为
\[
h_X=\frac{F_X'}{F_X},\qquad h_T=\frac{F_T'}{F_T}.
\]
条件分布
\[
\pi_t(dz)=\frac{F_X(t-z)f_Z(z)}{F_T(t)}\,dz,\qquad 0<z<t
\]
及
\[
q_t(z)=h_X(t-z),\qquad A_t(z)=-h_X'(t-z)
\]
给出精确恒等式
\[
h_T'=-\mathbb E_{\pi_t}A_t+\operatorname{Var}_{\pi_t}(q_t).
\]

本轮不再要求条件势 \(V_t=-\log F_X(t-z)-\log f_Z(z)\) 全局凸，而研究不依赖对数凹性的条件方差控制。

## 2. 真正需要的方差不等式

由于
\[
h_X'(t)=-A_t(0),
\]
有
\[
\boxed{
h_T'-h_X'
=A_t(0)-\mathbb E_{\pi_t}A_t
+\operatorname{Var}_{\pi_t}(q_t).
}
\]
因此，要证明相对 CDF 曲率
\[
(\log(F_T/F_X))''=h_T'-h_X'<0,
\]
真正需要的是
\[
\boxed{
\operatorname{Var}_{\pi_t}(q_t)
<
\mathbb E_{\pi_t}A_t-A_t(0).
}
\tag{W}
\]
单纯的 \(\operatorname{Var}(q_t)\le\mathbb E A_t\) 只重获 \(h_T'\le0\)，不能给出这个相对比较。

## 3. Stein kernel 与一维加权 Poincaré 表示

对条件密度 \(p_t\)，令
\[
\bar z_t=\int_0^t z p_t(z)\,dz.
\]
定义一维 Stein kernel
\[
\boxed{
\tau_t(z)=
\frac{1}{p_t(z)}
\int_z^t(y-\bar z_t)p_t(y)\,dy.
}
\]
等价地，也可使用下端积分的负号形式。对 \(0<z<t\)，\(\tau_t(z)\ge0\)。

对紧支撑区间上的绝对连续概率密度及相应 Sobolev 函数，一维 Stein-kernel 加权 Poincaré 不等式给出
\[
\boxed{
\operatorname{Var}_{\pi_t}(g)
\le
\mathbb E_{\pi_t}\!\left[\tau_t(Z)g'(Z)^2\right].
}
\]
取 \(g=q_t\)，利用 \(q_t'=A_t\)，得到
\[
\operatorname{Var}_{\pi_t}(q_t)
\le
\mathbb E_{\pi_t}[\tau_t A_t^2].
\]
故一个足以推出 (W) 的模型专属条件是
\[
\boxed{
\mathbb E_{\pi_t}[\tau_t(Z)A_t(Z)^2]
<
\mathbb E_{\pi_t}A_t-A_t(0).
}
\tag{S}
\]
此路线不要求 \(V_t''>0\)。条件密度 \(p_t\) 在 \((0,t)\) 为正且支撑有界，故其二阶矩有限；端点的密度消失可按弱导数/Sobolev 形式处理。

还有一个精确的双积分恒等式。令
\[
P_t(x)=\int_0^x p_t(z)\,dz,\qquad \bar P_t=1-P_t.
\]
则
\[
\boxed{
\begin{aligned}
\operatorname{Var}_{\pi_t}(q_t)
={}&\int_0^t\!\int_0^t A_t(x)A_t(y)\\
&\times P_t(x\wedge y)\bar P_t(x\vee y)\,dx\,dy.
\end{aligned}
}
\]
因此 (W) 也可写成明确的正核比较：
\[
\begin{aligned}
&\int_0^t\!\int_0^t A_t(x)A_t(y)
P_t(x\wedge y)\bar P_t(x\vee y)\,dx\,dy\\
&\qquad<
\int_0^t[A_t(z)-A_t(0)]p_t(z)\,dz.
\end{aligned}
\tag{K}
\]

Germain–Swan 的 Chen–Wang 变分公式还允许用选定的正权重估计最优加权 Poincaré 常数。对本问题，有限时间区间的后续工作可比较 (S)、(K) 以及相应的加权 Poincaré 常数，而不是恢复条件势凸性。

## 4. Keilson、Assaf 等人的定理条件及其适用边界

### 4.1 Keilson (1971)

J. Keilson, “Log-concavity and log-convexity in passage time densities of diffusion and birth-death processes,” *Journal of Applied Probability* **8** (1971), 391–398, [DOI](https://doi.org/10.2307/3211909)。

可公开检索的 Cambridge 记录提供摘要：文章讨论扩散过程与生灭过程中的若干关键 passage-time densities 和 sojourn-time densities 的对数凹/对数凸及相应单峰性质。当前可访问记录没有展示全文定理陈述，因此不把摘要扩写成比原文可核实内容更强的普遍定理。

### 4.2 Keilson (1981)：可核实的具体类别

J. Keilson, “On the unimodality of passage time densities in birth-death processes,” *Statistica Neerlandica* **35** (1981), 49–55, [DOI](https://doi.org/10.1111/j.1467-9574.1981.tb00710.x)。

Wiley 摘要明确给出：

1. 对任意 ergodic birth-death process，从反射态 0 到任意 level \(n\) 的 passage-time density \(T_{0n}\) 对数凹，因而强单峰。
2. 相邻状态之间的 \(T_{n,n+1}\) 和 \(T_{n+1,n}\) 密度完全单调，因而单峰。
3. 对任意两状态之间的 \(T_{mn}\)，证明其密度单峰；并且该密度可写为一个完全单调密度与一个 PF 密度的卷积，后者对数凹。

这些是**单个 passage-time density** 的形状结论。它们没有直接给出两个不同吸收端的 cause-specific cumulative capture functions 之比的单峰性。

### 4.3 Assaf–Shaked–Shanthikumar (1985)

D. Assaf, M. Shaked and J. G. Shanthikumar, “First-passage times with \(PF_r\) densities,” *Journal of Applied Probability* **22** (1985), 185–196, [DOI](https://doi.org/10.2307/3213758)。

摘要给出的充分条件是：

- 连续时间 Markov 链为有限状态；
- 可以进行 uniformization；
- uniformization 后的嵌入离散时间转移矩阵是 \(TP_r\)。

在这些条件下，从状态 0 到其它状态的 first-passage-time density 属于 \(PF_r\)。特别地，\(r=2\) 对应密度对数凹。该结果依赖有限状态 uniformization 和嵌入矩阵的总正性条件，不能直接当作连续空间分层 Robin 扩散的定理。

### 4.4 Kijima (1987) 的谱推广

M. Kijima, “Spectral structure of the first-passage-time densities for classes of Markov chains,” *Journal of Applied Probability* **24** (1987), 631–643, [DOI](https://doi.org/10.2307/3214095)。

摘要指出，其扩展对象包括时间可逆 Markov 链、skip-free Markov 链以及带吸收态的 birth-death processes，并利用生成元的谱表示得到相应首次通过时间密度。这是与本研究谱方法最接近的背景之一，但结论仍然针对单个 first-passage-time density。

### 4.5 对本模型的逐项区分

本模型的 \(A_L(t),A_R(t)\) 是两侧边界的累计捕获概率。它们的导数是两个 cause-specific capture subdensities；总吸收时间密度是二者之和，而目标量是
\[
\Delta(t)=A_L(t)-A_R(t).
\]
我们的证明对象还包括
\[
Q(t)=\frac{G_R(t)}{G_L(t)},\qquad
r(t)=\frac{G_R'(t)}{G_L'(t)}.
\]

因此：

- Keilson 的单个 \(T_{mn}\) 密度单峰，不等价于 \(A_L-A_R\) 的零点数上限。
- Assaf 等人的 \(PF_r\) 密度结论，不控制两个 competing absorbing outcomes 的累计概率差。
- Kijima 的谱表示可借鉴来分析单个吸收时间的谱结构，但不直接给出两个吸收端响应之比的曲率。
- 即使对模型构造 birth-death 离散近似，也仍需单独证明相应 \(TP_r\) 性质在极限中保持，并保留左右吸收端的分解；这不是现成定理自动提供的。

故这些文献提供重要的概率工具和方法背景，但尚未覆盖本研究的最终命题
\[
N_{\rm cross}(\Delta)\le2.
\]

## 5. Robin 谱给出的近打平晚时加权方差控制

令
\[
\eta_i=\lambda_{i+1}-c,\qquad
\rho_i=R_i-c,\qquad
\sigma_i=L_i-c.
\]
记
\[
\eta=\eta_2=\lambda_3-c.
\]
严格谱交错给出
\[
0<\eta<\rho_2\le\sigma_2,
\]
并且 \(\eta_i>\eta\) 对 \(i\ge3\)。

指数倾斜后的 Laplace 变换为
\[
\widehat f_{X_c}(s)
=\prod_{i\ge2}\frac{1+s/\rho_i}{1+s/\eta_i},
\]
\[
\widehat f_{T_c}(s)
=\prod_{i\ge2}\frac{1+s/\sigma_i}{1+s/\eta_i}.
\]

在简单极点、无首极点消去的非退化情形，最慢极点 \(-\eta\) 的留数为
\[
B_X=
\eta\left(1-\frac{\eta}{\rho_2}\right)
\prod_{i\ge3}
\frac{1-\eta/\rho_i}{1-\eta/\eta_i}>0,
\]
\[
B_T=
\eta\left(1-\frac{\eta}{\sigma_2}\right)
\prod_{i\ge3}
\frac{1-\eta/\sigma_i}{1-\eta/\eta_i}>0.
\]
两者之比具有一个关键的精确表达式：
\[
\boxed{
\frac{B_T}{B_X}
=
\prod_{i\ge2}
\frac{1-\eta/\sigma_i}{1-\eta/\rho_i}
=
\mathbb E[e^{\eta Z_c}]
=:\mathcal M_Z(\eta)>1.
}
\]
最后一个不等式来自 \(Z_c\) 非退化且其指数矩在 \(\eta<\rho_2\) 处有限。

由谱反演，
\[
h_X'(t)\sim-\eta B_Xe^{-\eta t},
\qquad
h_T'(t)\sim-\eta B_Te^{-\eta t}.
\]
于是
\[
\boxed{
h_T'(t)-h_X'(t)
\sim
-\eta B_X[\mathcal M_Z(\eta)-1]e^{-\eta t}<0.
}
\]

结合第 2 节的精确矩恒等式，得到
\[
\boxed{
\mathbb E_{\pi_t}A_t-A_t(0)
-\operatorname{Var}_{\pi_t}(q_t)
\sim
\eta B_X[\mathcal M_Z(\eta)-1]e^{-\eta t}>0.
}
\tag{L}
\]
因此，对每个满足上述简单极点和严格谱间隙条件的固定 Robin 参数，目标方差不等式 (W) 在充分晚时成立。

### 5.1 首模逆序趋于打平

令
\[
\varepsilon=c-a\downarrow0.
\]
阈值 \(t_*\) 满足 \(h_T(t_*)=\varepsilon\)，所以
\[
t_*=
\frac1\eta\log\frac{B_T}{\varepsilon}+o(1).
\]
对任意固定 \(d>0\)，取
\[
t=t_*-d.
\]
则
\[
e^{-\eta t}\sim
\frac{\varepsilon e^{\eta d}}{B_T}.
\]
代入 (L)，得到
\[
\boxed{
\begin{aligned}
&\mathbb E_{\pi_t}A_t-A_t(0)
-\operatorname{Var}_{\pi_t}(q_t)\\
&\qquad\sim
\eta\varepsilon e^{\eta d}
\left(1-\frac1{\mathcal M_Z(\eta)}\right)>0.
\end{aligned}
}
\]
这说明在近打平的晚时尺度上，平均曲率差具有 \(O(\varepsilon)\) 的正主项，且严格压过条件方差。点态条件势凸性可以失败，而平均条件方差不等式仍然成立。

若需给出对 \(\kappa\uparrow\kappa_0\) 一致的版本，需要使用临界几何附近的统一谱隙和留数界；严格交错保证主极点不与第二极点合并。

## 6. 本轮结论与剩余工作

已经得到：

1. Keilson、Assaf 等人的相关定理均针对单个 passage-time density；其条件和结论与两个吸收端累计概率差的反超上限不同。
2. Stein kernel 给出不要求 \(V_t\) 凸的加权方差估计。
3. Robin 谱留数给出 \(B_T/B_X=\mathbb E e^{\eta Z_c}>1\) 的精确关系。
4. 在充分晚时，已经证明目标方差不等式 (W) 成立；在首模趋于打平的 \(t=t_*-d\) 尺度上，其正裕量有明确的 \(O(\varepsilon)\) 渐近。

尚未完成的是中间时间区域的统一估计，以及把 \(h_T'-h_X'\) 的比较扩展到完整
\[
\mathcal C(t)=(\log r(t))''.
\]
因此本轮结果不宣称 \(r\) 已经单峰，也不宣称最终 \(N_{\rm cross}(\Delta)\le2\) 已经证明。

## 7. 参考文献

1. Keilson, J. (1971). Log-concavity and log-convexity in passage time densities of diffusion and birth-death processes. *J. Appl. Probab.* 8, 391–398. [DOI](https://doi.org/10.2307/3211909).
2. Keilson, J. (1981). On the unimodality of passage time densities in birth-death processes. *Statistica Neerlandica* 35, 49–55. [DOI](https://doi.org/10.1111/j.1467-9574.1981.tb00710.x).
3. Assaf, D., Shaked, M. and Shanthikumar, J. G. (1985). First-passage times with \(PF_r\) densities. *J. Appl. Probab.* 22, 185–196. [DOI](https://doi.org/10.2307/3213758).
4. Kijima, M. (1987). Spectral structure of the first-passage-time densities for classes of Markov chains. *J. Appl. Probab.* 24, 631–643. [DOI](https://doi.org/10.2307/3214095).
5. Germain, G. and Swan, Y. (2023). A note on one-dimensional Poincaré inequalities by Stein-type integration. *Bernoulli* 29, 1714–1740. [DOI](https://doi.org/10.3150/22-BEJ1518).
