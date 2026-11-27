# 条件势点态凸性条件的谱参数障碍

## 1. 结论

当前拟证明的条件
\[
-h_X'(t-z)>(\log f_{Z_c}(z))'',
\qquad 0<z<t<t_*,
\]
等价于条件势
\[
V_t(z)=-\log F_X(t-z)-\log f_{Z_c}(z)
\]
在整个条件区间严格凸。

本轮发现：**这个点态凸性条件不能作为整个首模唯一逆序参数域上的统一引理。** 在一个具体的三层 Fick–Robin 几何族中，可以令 \(c-a\downarrow0\)，使 \(t_*\to\infty\)；与此同时，\(f_{Z_c}\) 在一个固定的晚时间点具有严格正的对数曲率，而 \(-h_X'(t-z)\) 可以趋于零。因此点态不等式最终失败。

这否定的是 Brascamp–Lieb 全域证明路线所需的充分条件，不是否定 \(r(t)\) 单峰性，也不构成累计捕获反超次数上限的反例。

## 2. 一个物理 Robin 参数族

保持原有三层几何：
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区间长度 \(1/2\)、扩散率 \(D_R=1\)，两端共同 Robin 参数为 \(\kappa\)。

右半第一 Robin–Dirichlet 波数 \(k\in(\pi,2\pi)\) 满足
\[
\cos\frac{k}{2}+\frac{\kappa}{k}\sin\frac{k}{2}=0,
\qquad
\kappa=-k\cot\frac{k}{2}.
\]
左半特征函数由两层传递矩阵
\[
M_D(w;k)=
\begin{pmatrix}
\cos(kw/\sqrt D)&\dfrac{\sin(kw/\sqrt D)}{\sqrt D\,k}\\[5pt]
-\sqrt D\,k\sin(kw/\sqrt D)&\cos(kw/\sqrt D)
\end{pmatrix}
\]
给出。令 \(M=M_{D_2}(w_2;k)M_{D_1}(w_1;k)\)，左半 Robin–Dirichlet 特征方程为
\[
Y_L(k;\kappa)=M_{11}(k)+\kappa M_{12}(k)=0.
\]

定义 \(\kappa_0\) 使左右第一半谱相等：
\[
L_1(\kappa_0)=R_1(\kappa_0)=k_0^2.
\]
等价地，\(k_0\) 是
\[
Y_L\!\left(k; -k\cot\frac{k}{2}\right)=0
\]
在 \((\pi,2\pi)\) 内的根。

### 区间认证

60 位 mpmath 区间算术给出
\[
k_0\in
[5.3476625759351,\ 5.3476625759352],
\]
\[
\boxed{\kappa_0\in
[10.5862272293954682,\ 10.5862272293969815]},
\]
以及
\[
L_1=R_1=k_0^2
\in[28.59749502605682,\ 28.59749502605790].
\]

交叉特征函数在波数区间两端严格异号，其导数区间为
\[
\frac{d}{dk}Y_L\!\left(k;-k\cot\frac{k}{2}\right)
\in[3.3486808264382,\ 3.3486808264404]>0,
\]
且 \(d\kappa/dk>0\)。因此在 \(\kappa_0\) 左侧充分小的邻域，\(L_1<R_1\)。由全谱与半谱的 Sturm 交错，
\[
a=L_1<c=\lambda_2<b=R_1,
\qquad
\varepsilon:=c-a>0,
\]
而
\[
\boxed{\varepsilon\longrightarrow0
\quad(\kappa\uparrow\kappa_0).}
\]

在 \(\kappa=\kappa_0\)，左右半区第一模态都在中点为零，可按通量匹配拼接为全区间具有一个内节点的特征函数。因此其特征值正是 \(\lambda_2\)，即 \(c=a=b\)。

### 高阶谱排序

在 \(\kappa_0\) 的区间根证书为：

| 谱根 | 波数平方区间 |
|---|---:|
| \(R_2\) | [120.13933470393, 120.13933470459] |
| \(L_2\) | [177.40882851926, 177.40882851980] |
| \(R_3\) | [283.27822208445, 283.27822208580] |
| \(L_3\) | [665.12864222663, 665.12864222767] |
| \(R_4\) | [522.51412236459, 522.51412236598] |
| \(L_4\) | [1424.05281670960, 1424.05281671112] |

特别地，
\[
R_2<L_2<R_3<L_3.
\]

还可以解析控制所有高模。左半光学长度和右半光学长度分别为
\[
\ell_L=\frac{w_1}{5}+2w_2=\frac7{22},
\qquad
\ell_R=\frac12,
\qquad
\delta_\ell=\ell_R-\ell_L=\frac2{11}.
\]
左半 Prüfer 相位可写为
\[
\Theta_L(k)=k\ell_L-\arctan\frac{\kappa}{5k}+\delta_L(k),
\qquad |\delta_L(k)|<\frac{\pi}{2},
\]
其第 \(n\) 个根满足 \(\Theta_L(k_{L,n})=(n-\tfrac12)\pi\)。右半相位为
\[
\Theta_R(k)=\frac{k}{2}+\arctan\frac{k}{\kappa},
\qquad \Theta_R(k_{R,n})=n\pi.
\]
因此
\[
\Theta_L(k_{R,n})+\frac{\pi}{2}-n\pi
<
-k_{R,n}\delta_\ell+\pi.
\]
又 \(k_{R,n}>(2n-1)\pi\)。当 \(n\ge4\) 时，
\[
k_{R,n}\delta_\ell>7\pi\frac2{11}
=\frac{14\pi}{11}>\pi,
\]
故 \(\Theta_L(k_{R,n})+\pi/2<n\pi\)。由 Prüfer 相位严格递增，
\[
L_n>R_n,\qquad n\ge4.
\]
结合 \(n=2,3\) 的区间根比较，得到
\[
\boxed{L_i>R_i\quad\text{对所有 }i\ge2.}
\]
因此 \(\kappa<\kappa_0\) 且充分接近 \(\kappa_0\) 时，这一族确实处于首模唯一逆序、高阶谱正常排序的参数域。

## 3. \(Z_c\) 密度具有固定的晚时正曲率区

记
\[
\rho_i=R_i-c,\qquad \sigma_i=L_i-c.
\]
指数倾斜只给密度乘上指数因子：
\[
f_{Z_c}(z)=\frac{e^{cz}}{M_Z}f_\nu(z).
\]
所以
\[
\boxed{(\log f_{Z_c}(z))''=(\log f_\nu(z))''.}
\]
其对数曲率不依赖倾斜参数 \(c\)。

对谱乘积
\[
\widehat f_{Z_c}(s)
=\prod_{i\ge2}\frac{1+s/\sigma_i}{1+s/\rho_i},
\]
极点 \(-\rho_j\) 的留数为
\[
A_j=\rho_j\left(1-\frac{\rho_j}{\sigma_j}\right)
\prod_{\substack{i\ge2\\i\ne j}}
\frac{1-\rho_j/\sigma_i}{1-\rho_j/\rho_i}.
\]

在本参数族的极限点，
\[
\rho_2<\sigma_2<\rho_3,\qquad \sigma_3>\rho_3,
\]
且 \(L_i>R_i\) 对 \(i\ge4\)。因此 \(A_2,A_3>0\)。正光学长度差保证 Laplace 反演可在前三个极点之后移动，得到
\[
f_{Z_c}(z)
=A_2e^{-\rho_2z}+A_3e^{-\rho_3z}
+O(e^{-\eta z}),
\qquad \rho_3<\eta<\rho_4.
\]
从而
\[
\boxed{
(\log f_{Z_c}(z))''
\sim
\frac{A_3}{A_2}(\rho_3-\rho_2)^2
e^{-(\rho_3-\rho_2)z}>0.
}
\]
所以存在一个固定的有限 \(z_0>0\)，使极限几何满足
\[
g_0:=(\log f_{Z_{c_0}}(z_0))''>0.
\]
谱根的连续依赖及正时间热核的局部一致收敛保证，\(\kappa<\kappa_0\) 充分接近时仍有
\[
(\log f_{Z_c}(z_0))''>\frac{g_0}{2}>0.
\]

## 4. \(h_X'\) 的晚时衰减与 \(t_*\) 发散

记
\[
\gamma_i=\lambda_{i+1}-c.
\]
由谱交错，
\[
0<\gamma_2<\rho_2\le\sigma_2,
\qquad
\gamma_i>\gamma_2\quad(i\ge3).
\]

\(X_c\) 和 \(T_c=X_c+Z_c\) 的 Laplace 变换分别为
\[
K(s)=\prod_{i\ge2}\frac{1+s/\rho_i}{1+s/\gamma_i},
\qquad
K(s)\widehat f_{Z_c}(s)
=\prod_{i\ge2}\frac{1+s/\sigma_i}{1+s/\gamma_i}.
\]
它们在最慢极点 \(-\gamma_2\) 的留数分别为
\[
B_X=\gamma_2\left(1-\frac{\gamma_2}{\rho_2}\right)
\prod_{i\ge3}
\frac{1-\gamma_2/\rho_i}{1-\gamma_2/\gamma_i}>0,
\]
\[
B_T=\gamma_2\left(1-\frac{\gamma_2}{\sigma_2}\right)
\prod_{i\ge3}
\frac{1-\gamma_2/\sigma_i}{1-\gamma_2/\gamma_i}>0.
\]
故
\[
-h_X'(u)\sim\gamma_2B_Xe^{-\gamma_2u},
\qquad
h_T(t)\sim B_Te^{-\gamma_2t}.
\]

令 \(\varepsilon=c-a\)，阈值时间满足 \(h_T(t_*)=\varepsilon\)。于是
\[
\boxed{
t_*=\frac1{\gamma_2}\log\frac{B_T}{\varepsilon}+o(1)
\longrightarrow\infty
\quad(\kappa\uparrow\kappa_0).
}
\]
这里 \(\gamma_2,B_X,B_T\) 在极限点均有严格正的有限极限。

## 5. 点态曲率补偿条件失败

固定上面选取的 \(z_0\)，再固定任意 \(d>0\)，取
\[
t_\kappa=t_*(\kappa)-d.
\]
当 \(\kappa\uparrow\kappa_0\) 时，
\[
0<z_0<t_\kappa<t_*(\kappa),
\qquad
t_\kappa-z_0\longrightarrow\infty.
\]
利用晚时渐近和 \(e^{-\gamma_2t_*}\sim\varepsilon/B_T\)，
\[
\begin{aligned}
-h_X'(t_\kappa-z_0)
&\sim \gamma_2B_X
e^{-\gamma_2(t_*-d-z_0)}\\
&\sim
\frac{\gamma_2B_X}{B_T}
e^{\gamma_2(d+z_0)}\varepsilon
\longrightarrow0.
\end{aligned}
\]
但
\[
(\log f_{Z_c}(z_0))''>\frac{g_0}{2}>0.
\]
所以对所有充分接近 \(\kappa_0\) 的首模唯一逆序参数，
\[
\boxed{
-h_X'(t_\kappa-z_0)
<
(\log f_{Z_c}(z_0))''.
}
\]

即拟证明的点态条件在其指定区域内确实失败。相应地，
\[
V_{t_\kappa}''(z_0)<0.
\]

## 6. 对原研究目标的影响

1. 全局条件势严格凸性不能作为整个首模唯一逆序参数域的统一引理。
2. Brascamp–Lieb 可以继续用于 \((\log f_{Z_c})''\le0\) 的早时区域，或用于点态补偿确实成立的受限参数区。
3. 全域反超次数定理仍然开放；本结果不推出 \(r(t)\) 非单峰，也不构成 \(N_{\rm cross}(\Delta)>2\) 的物理反例。
4. 下一步应改用条件方差的加权 Poincaré/Stein-kernel 控制，或直接研究 \(N/D\) 的谱卷积比值，而不再要求整个条件势逐点凸。

## 7. 相关前人工作

- Sumita and Masuda, “Classes of Probability Density Functions Having Laplace Transforms with Negative Zeros and Poles,” *Advances in Applied Probability* (1987), [DOI](https://doi.org/10.2307/1427410)。负实极点/零点及指数卷积密度的分类背景。
- Keilson, “Log-concavity and Log-convexity in Passage Time Densities of Diffusion and Birth-Death Processes,” *Journal of Applied Probability* (1971), [DOI](https://doi.org/10.2307/3211909)。首次通过时间密度的对数曲率背景。
- Saumard and Wellner, “Log-Concavity and Strong Log-Concavity: A Review,” *Statistics Surveys* (2014), [DOI](https://doi.org/10.1214/14-SS107)。Prékopa 与 Brascamp–Lieb 等工具的标准参考。

本轮的新增内容是：将谱留数的晚时对数凸性与 \(c-a\downarrow0\) 的物理 Robin 参数族结合，证明了全域点态曲率补偿条件失效。并未把这一障碍误写成最终捕获反超定理的反例。

## 8. 计算与归档记录

- 精确输入：curvature_compensation_obstruction_input.json
- 区间认证源码：curvature_compensation_obstruction.py
- 原生区间输出：curvature_compensation_obstruction_output.txt
- 独立高精度复核源码：curvature_compensation_independent_check.py
- 独立复核输出：curvature_compensation_independent_output.txt
- 区间运算：Python 3.13.5，mpmath 1.3.0，60 位十进制区间精度。
- 独立复核：70 位 mpmath 高精度根求解。
- 环境：Linux 6.18.57 x86_64，AMD EPYC 9V74；未使用 GPU。
- 输入 SHA-256：1e513a24125aa906d51b8087dbc030763aeb496939db3a3961af53e53d4fcbe2
- 区间源码 SHA-256：7c832226e689691e0039d35d7ba1169b4464c75cbbbf6eacf30406571f88e816
- 区间输出 SHA-256：74d333499d47736ac00a8ad2c5f63d223217fd8240981edf3278a57d6141868b
- 独立源码 SHA-256：5bda669a3d0ceb19d3aff3c47f109db4a30f7519907a62629160b3d8b6047389
- 独立输出 SHA-256：e81b65385ef37bd58e144ddb9abfe99cee1aea1f955580a414dfbd1089c70d6f

区间证书的根括区端点严格异号，交叉导数和 \(d\kappa/dk\) 的区间均严格为正；独立高精度求根与区间结果一致。高阶谱排序由上述 Prüfer 相位界与有限根区间共同控制。
