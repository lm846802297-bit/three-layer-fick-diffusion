# Robin score 曲率的高阶谱项聚合控制

## 1. 目标

考虑静态打平 Robin 基准：
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区长度 \(1/2\)、扩散率 \(D_R=1\)，共同 Robin 参数 \(\kappa>0\)，并令
\[
c=L_1=R_1=\lambda_2.
\]
定义
\[
\rho_j=R_j-c,\qquad \sigma_j=L_j-c,\qquad j\ge2.
\]
延迟密度的 Laplace 变换为
\[
\Phi_Z(s)=\prod_{j\ge2}\frac{1+s/\sigma_j}{1+s/\rho_j}.
\]
在留数展开成立时，
\[
f_Z(t)=\sum_{j\ge2}A_je^{-\rho_jt},
\]
且
\[
Q_Z(t):=f_Zf_Z''-(f_Z')^2
=\sum_{2\le i<j}C_{ij}e^{-\lambda_{ij}t},
\]
其中
\[
C_{ij}=A_iA_j(\rho_j-\rho_i)^2,\qquad
\lambda_{ij}=\rho_i+\rho_j.
\]
因为 \(f_Z>0\)，\(\operatorname{sgn}(\log f_Z)''=\operatorname{sgn}Q_Z\)。

本节建立一个模型专属的高阶留数增长界，并证明双谱尾在任何固定正时间区间上可由有限项以显式 Gaussian 误差统一控制。它不要求高阶留数同号。

## 2. 两个半谱特征函数

令 \(k=\sqrt\lambda\)。按左半区从外壁向中点传播，取初值 \(y(0)=1\)、\(D y'(0)=\kappa\)，左半区到中点的特征函数为
\[
\begin{aligned}
Y_L(k)={}&
\cos\frac{5k}{66}\cos\frac{16k}{66}
-10\sin\frac{5k}{66}\sin\frac{16k}{66}\\
&+\frac{\kappa}{k}\left[
\frac15\cos\frac{16k}{66}\sin\frac{5k}{66}
+2\sin\frac{16k}{66}\cos\frac{5k}{66}
\right].
\end{aligned}
\tag{2.1}
\]
右半区是均匀介质，其特征函数为
\[
Y_R(k)=\cos\frac{k}{2}+\frac{\kappa}{k}\sin\frac{k}{2}.
\tag{2.2}
\]
令 \(F_L(\lambda)=Y_L(\sqrt\lambda)/Y_L(0)\)、\(F_R(\lambda)=Y_R(\sqrt\lambda)/Y_R(0)\)。本基准中
\[
Y_L(0)=Y_R(0)=1+\frac{\kappa}{2}.
\]
\(F_L,F_R\) 是相应半谱的归一化规范乘积。它们在 \(\lambda=c\) 有共同简单零点。约去共同首因子后，
\[
\Phi_Z(s)=C_0
\frac{\bar F_L(c-s)}{\bar F_R(c-s)},
\qquad
C_0=\frac{F_R'(c)}{F_L'(c)}>0.
\tag{2.3}
\]
这里撇号表示对 \(\lambda\) 求导，\(\bar F_{L,R}(\lambda)=F_{L,R}(\lambda)/(1-\lambda/c)\)。

在右半第 \(j\) 个根 \(k_j=\sqrt{R_j}\) 处，留数因而有精确表达式
\[
\boxed{
A_j=-2C_0k_j\frac{Y_L(k_j)}{Y_R'(k_j)},
}
\tag{2.4}
\]
其中 \(Y_R'\) 表示对 \(k\) 求导。

## 3. 留数只有线性增长

右半根方程为
\[
\cos(k_j/2)+\frac{\kappa}{k_j}\sin(k_j/2)=0.
\]
因此
\[
(2j-1)\pi<k_j<2j\pi.
\tag{3.1}
\]
在根处代入 \(\cos(k_j/2)=-(\kappa/k_j)\sin(k_j/2)\)，得到
\[
Y_R'(k_j)=
-\sin\frac{k_j}{2}
\left[
\frac12+\frac{\kappa^2/2+\kappa}{k_j^2}
\right].
\]
又有
\[
\left|\sin\frac{k_j}{2}\right|
=\frac{k_j}{\sqrt{k_j^2+\kappa^2}},
\]
所以
\[
\boxed{
|Y_R'(k_j)|
\ge\frac{k_j}{2\sqrt{k_j^2+\kappa^2}}.
}
\tag{3.2}
\]

由 (2.1) 的传递矩阵表达式，
\[
|Y_L(k)|
\le 1+10+\frac{\kappa}{k}\left(\frac15+2\right)
=11+\frac{11\kappa}{5k}.
\]
令 \(k_2=\sqrt{R_2}\)，则对所有 \(j\ge2\)，
\[
|Y_L(k_j)|\le M_L,\qquad
M_L:=11+\frac{11\kappa}{5k_2}.
\]
结合 (2.4)、(3.2)，
\[
|A_j|
\le
4C_0M_L\sqrt{1+\frac{\kappa^2}{k_2^2}}\,k_j.
\]
再用 \(k_j<2\pi j\)，得到

\[
\boxed{
|A_j|\le B j,\qquad
B:=8\pi C_0M_L
\sqrt{1+\frac{\kappa^2}{k_2^2}},
\quad j\ge2.
}
\tag{3.3}
\]

因此，尽管 \(A_j\) 的符号无限次改变，其绝对值只有线性增长。

由 (3.1) 及本基准 \(c<5\pi^2\)，还得到
\[
\begin{aligned}
\rho_j
&=k_j^2-c\\
&>(2j-1)^2\pi^2-c\\
&\ge \pi^2j^2,
\end{aligned}
\]
以及
\[
\rho_j<4\pi^2j^2.
\]
故
\[
\boxed{
\pi^2j^2\le\rho_j\le4\pi^2j^2,\qquad j\ge2.
}
\tag{3.4}
\]

## 4. 双谱高阶项的显式聚合界

定义
\[
Q_N(t)=
\sum_{2\le i<j\le N}
C_{ij}e^{-\lambda_{ij}t},
\]
以及正项和
\[
S_m(t)=\sum_{j\ge2}|A_j|\rho_j^m e^{-\rho_jt},
\qquad
S_{m,N}(t)=\sum_{j>N}|A_j|\rho_j^m e^{-\rho_jt}.
\]
当 \(t>0\) 时，这些级数均绝对收敛。

由
\[
(\rho_j-\rho_i)^2\le2(\rho_i^2+\rho_j^2)
\]
和 \(i<j\)、\(j>N\)，得到
\[
\boxed{
|Q_Z(t)-Q_N(t)|
\le
2\left[S_0(t)S_{2,N}(t)+S_2(t)S_{0,N}(t)\right].
}
\tag{4.1}
\]

对 \(Q_Z'\) 求导，并使用
\[
(\rho_j-\rho_i)^2(\rho_i+\rho_j)
\le2(\rho_i^3+\rho_i^2\rho_j+\rho_i\rho_j^2+\rho_j^3),
\]
得到
\[
\boxed{
\begin{aligned}
|Q_Z'(t)-Q_N'(t)|
\le2\big[&
S_3(t)S_{0,N}(t)+S_2(t)S_{1,N}(t)\\
&+S_1(t)S_{2,N}(t)+S_0(t)S_{3,N}(t)
\big].
\end{aligned}
}
\tag{4.2}
\]

由 (3.3)–(3.4)，对任意固定 \(\tau>0\) 及 \(t\ge\tau\)，
\[
S_{m,N}(t)\le U_{m,N}(\tau),
\]
其中完全显式的上界为
\[
\boxed{
U_{m,N}(\tau)
=
B(4\pi^2)^m
\sum_{j>N}j^{2m+1}e^{-\pi^2j^2\tau}.
}
\tag{4.3}
\]
类似地，令 \(U_m(\tau)\) 为将求和下限改成 \(j=2\) 的总和上界。于是
\[
\boxed{
\sup_{t\ge\tau}|Q_Z-Q_N|
\le
2[U_0U_{2,N}+U_2U_{0,N}],
}
\tag{4.4}
\]
以及
\[
\boxed{
\begin{aligned}
\sup_{t\ge\tau}|Q_Z'-Q_N'|
\le2[&
U_3U_{0,N}+U_2U_{1,N}\\
&+U_1U_{2,N}+U_0U_{3,N}].
\end{aligned}
}
\tag{4.5}
\]
所有 \(U\) 均在 \(\tau\) 处取值。

对固定 \(\tau>0\)，高阶尾满足
\[
U_{m,N}(\tau)
=O_\tau\!\left(N^{2m}e^{-\pi^2\tau N^2}\right).
\tag{4.6}
\]
所以
\[
\boxed{
Q_N\longrightarrow Q_Z
\quad\text{在 }[\tau,\infty)\text{ 上以显式误差在 }C^1\text{ 中一致收敛。}
}
\tag{4.7}
\]

这是本轮的主要聚合结果。它控制的是所有高阶成对项的绝对总贡献，不要求也不试图消除单个 \(A_j\) 的符号振荡。

## 5. 为什么小时间必须保留整个高频块

当 \(\tau\downarrow0\)，(4.3) 的尾界不再小。这不是估计缺陷，而是小时间渐近的真实结构：有效参与的谱模数随 \(t\downarrow0\) 增长，任何固定有限截断都不能一致描述原点附近的曲率。

这一点可以从同一个传递矩阵直接看出。令 \(q=\sqrt{s-c}\)。在 \(q\to+\infty\) 时，
\[
Y_L(iq)=
\frac{11}{4}e^{q\ell_L}
\left(1+\frac{\kappa}{5q}
+O(e^{-2\ell_{\min}q})\right),
\]
\[
Y_R(iq)=
\frac12e^{q\ell_R}
\left(1+\frac{\kappa}{q}
+O(e^{-2\ell_Rq})\right),
\]
其中
\[
\ell_L=\frac7{22},\qquad
\ell_R=\frac12,\qquad
\delta_\ell=\ell_R-\ell_L=\frac2{11},
\qquad
\ell_{\min}=\frac5{66}.
\]
代入 (2.3)，并将 \(q=\sqrt{s-c}\) 展开为 \(\sqrt s\)，得到
\[
\boxed{
\Phi_Z(s)=
\frac{11C_0}{2}e^{-\delta_\ell\sqrt s}
\left[
1+\frac{\gamma_1}{\sqrt s}+O(s^{-1})
\right],
}
\tag{5.1}
\]
其中
\[
\boxed{
\gamma_1=\frac{\delta_\ell c}{2}-\frac{4\kappa}{5}.
}
\tag{5.2}
\]
这里的指数型传递矩阵余项在逆 Laplace 鞍点下产生比主项更短的指数尺度。

利用
\[
\mathcal L^{-1}\{e^{-\delta\sqrt s}\}
=
\frac{\delta}{2\sqrt\pi}t^{-3/2}
e^{-\delta^2/(4t)}
\]
及
\[
\mathcal L^{-1}\{s^{-1/2}e^{-\delta\sqrt s}\}
=
\frac1{\sqrt{\pi t}}e^{-\delta^2/(4t)},
\]
得到带可微渐近余项的形式
\[
\boxed{
f_Z(t)=
\frac{11C_0\delta_\ell}{4\sqrt\pi}
t^{-3/2}e^{-\delta_\ell^2/(4t)}
\left[
1+\frac{2\gamma_1}{\delta_\ell}t+O(t^2)
\right].
}
\tag{5.3}
\]
因此
\[
\boxed{
(\log f_Z(t))''
=
-\frac{\delta_\ell^2}{2t^3}
+\frac{3}{2t^2}+O(1)<0
}
\tag{5.4}
\]
对充分小的正时间成立。

(5.1)–(5.4) 描述的是无限高阶谱项的集体贡献；它不能由固定 \(N\) 的 \(Q_N\) 在 \(t=0\) 附近代替。

## 6. 晚时间的首对谱项控制

已有
\[
A_2>0,\qquad A_3>0,\qquad \rho_2<\rho_3.
\]
因此首对系数
\[
C_{23}=A_2A_3(\rho_3-\rho_2)^2>0
\]
对应最小指数率
\[
\lambda_{23}=\rho_2+\rho_3.
\]
所有其它成对谱率均严格大于 \(\lambda_{23}\)。

定义
\[
\mathcal E(T)=
\sum_{\substack{i<j\\(i,j)\ne(2,3)}}
\frac{|A_iA_j|(\rho_j-\rho_i)^2}{C_{23}}
e^{-[\rho_i+\rho_j-\rho_2-\rho_3]T}.
\tag{6.1}
\]
由 (4.4) 所给的绝对收敛性，对任意 \(T>0\) 该和有限；又由所有指数差严格为正，
\[
\mathcal E(T)\longrightarrow0,\qquad T\to\infty.
\]
所以存在有限 \(T\)，使 \(\mathcal E(T)<1\)，并且
\[
\boxed{
Q_Z(t)>0,\qquad t\ge T.
}
\tag{6.2}
\]
这给出首对谱项支配的定量版本，不依赖高阶留数符号。

## 7. 低模核心与唯一零点的剩余任务

现有静态打平根区间给出
\[
R_2<L_2<R_3<R_4<L_3.
\]
由留数符号公式
\[
\operatorname{sgn}A_j=(-1)^{m(\rho_j+)-1}
\]
得到
\[
A_2>0,\qquad A_3>0,\qquad A_4<0.
\]
因此前四模的曲率分子为
\[
Q_4(t)=
C_{23}e^{-(\rho_2+\rho_3)t}
-C_{24}^{\!*}e^{-(\rho_2+\rho_4)t}
-C_{34}^{\!*}e^{-(\rho_3+\rho_4)t},
\]
其中 \(C_{24}^{\!*},C_{34}^{\!*}>0\)。乘以 \(e^{(\rho_2+\rho_3)t}\) 后，
\[
e^{(\rho_2+\rho_3)t}Q_4(t)
=
C_{23}
-C_{24}^{\!*}e^{-(\rho_4-\rho_3)t}
-C_{34}^{\!*}e^{-(\rho_4-\rho_2)t}
\]
严格递增。因此 \(Q_4\) 至多有一个零点。

但 \(Q_4\) 不能直接代替完整 \(Q_Z\)：在小时间，高阶谱块不可忽略；在中间时间，需使用 (4.4)–(4.5) 控制其聚合误差。

完整唯一零点证明现在被精确分为三个区间任务：

1. **小时间：** 用 (5.1) 的完整高频传递矩阵鞍点展开，建立一个显式 \(t_{\rm small}>0\)，使 \(Q_Z<0\) 于 \((0,t_{\rm small}]\)。
2. **中间时间：** 选取 \(N\)，用 (4.4)–(4.5) 控制高阶尾；对有限 \(Q_N\) 建立一个单零点区间，并验证根附近 \(Q_N'>0\)、其余部分有严格符号裕量。
3. **晚时间：** 用 (6.1) 选择 \(T\)，使 \(Q_Z>0\) 对全部 \(t\ge T\) 成立。

本轮已经完成高阶留数的统一增长界和双谱尾的显式 \(C^1\) 聚合界；尚未完成中间时间有限和的唯一零点认证。因此尚不宣称 \(G_Z\) 唯一变号。

## 8. 本轮结论

\[
\boxed{
\begin{gathered}
|A_j|\le Bj,\qquad \rho_j\asymp j^2,\\
Q_N\to Q_Z\quad\text{在每个 }[\tau,\infty),\ \tau>0,\text{ 上以显式 }C^1\text{ 误差收敛},\\
G_Z(t)<0\text{ 对充分小 }t>0,\qquad
G_Z(t)>0\text{ 对充分大 }t.
\end{gathered}
}
\]
无限次留数符号振荡现在有了一个可控的整体处理方式：小时间保留完整高频传递矩阵，固定正时间区间使用有限谱和加 Gaussian 尾界，晚时间使用首对极点支配。剩余的核心工作是中间有限谱和的唯一零点证明。
