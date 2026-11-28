# Robin score 曲率的高阶谱项聚合控制

## 1. 目标

考虑静态打平 Robin 基准：
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区长度 \(1/2\)、扩散率 \(D_R=1\)，共同 Robin 参数 \(\kappa>0\)，并令
\[
c=\lambda_2,\qquad L_1<c<R_1.
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
\(F_L,F_R\) 是相应半谱的归一化规范乘积。这里 \(c=\lambda_2\) 是全区间第二 Robin 本征值，并不是两个半谱的共同根；一般有 \(L_1<c<R_1\)。定义去掉各自首根后的半谱乘积
\[
\bar F_L(\lambda)=\frac{F_L(\lambda)}{1-\lambda/L_1},
\qquad
\bar F_R(\lambda)=\frac{F_R(\lambda)}{1-\lambda/R_1}.
\]
于是
\[
\Phi_Z(s)=C_0\frac{\bar F_L(c-s)}{\bar F_R(c-s)},
\qquad
C_0=\frac{\bar F_R(c)}{\bar F_L(c)}>0,
\tag{2.3}
\]
其中 \(C_0\) 由 \(\Phi_Z(0)=1\) 固定。

在右半第 \(j\) 个根 \(k_j=\sqrt{R_j}\) 处，留数因而有精确表达式
\[
\boxed{
A_j=-2C_0k_j\frac{Y_L(k_j)}{Y_R'(k_j)}
\frac{1-R_j/R_1}{1-R_j/L_1}.
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
对 \(j\ge2\)，有 \(R_j>R_1>L_1\)，从而
\[
0<\frac{1-R_j/R_1}{1-R_j/L_1}<1.
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
\frac{11C_0L_1}{2R_1}e^{-\delta_\ell\sqrt s}
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
\frac{11C_0L_1\delta_\ell}{4R_1\sqrt\pi}
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

结合后续归档的短时反射路径证书、Q30 短桥、Q50 外侧区间证书及 Q30 中间窗证书，当前 \\(\\kappa=4\\) 静态打平基准的时间轴已分段覆盖：

1. **短时解析区：** \\(0<t\\leq0.0015\\)，由显式 Bromwich 反射余项界证明 \\(Q_Z<0\\)。
2. **短桥：** \\(0.0015\\leq t\\leq0.004\\)，由 Q30、R2–R70 根括区及有限/无限谱尾界证明 \\(Q_Z<0\\)。
3. **中间负曲率区：** \\(0.004\\leq t\\leq0.010\\)，由 Q50 Taylor-box 证书证明 \\(Q_Z<0\\)。
4. **唯一穿越区：** \\(0.010\\leq t\\leq0.012\\)，由 Q30 的 C1 区间证书证明 \\(Q_Z'>0\\)，且唯一零点位于 \\(0.01102,0.01103\\) 之间。
5. **晚时正曲率区：** \\(t\\geq0.012\\)，由首对谱项支配证明 \\(Q_Z>0\\)。

因此，对这一具体基准，\\((\\log f_Z)''\\) 的全时间唯一简单零点结论已经闭合。该结论不自动推广到其它 \\(\\kappa\\) 或几何参数。

## 8. 当前全时间结论

结合短时反射路径证书、Q30 短桥、Q50 外侧区间证书、Q30 中间窗口证书和首对谱项支配，当前 \(\kappa=4\) 静态打平基准的全时间符号已经覆盖：
\[
\boxed{
Q_Z(t)
\begin{cases}
<0,&0<t<z_*,\\
=0,&t=z_*,\\
>0,&t>z_*,
\end{cases}
\qquad z_*\in(0.01102,0.01103).
}
\]
因为 \(f_Z(t)>0\)，这等价于 \((\log f_Z)''\) 恰有一个简单零点，并由负变正。

具体证书分区为：
\[
\begin{aligned}
(0,0.0015]&:\quad\text{短时反射路径解析界},\\
[0.0015,0.010]&:\quad\text{Q50 Taylor-box 与 Gaussian 双谱尾},\\
[0.010,0.012]&:\quad\text{Q30 的 }C^1\text{ 区间单穿越证书},\\
[0.012,\infty)&:\quad\text{首对谱项严格支配}.
\end{aligned}
\]
Q30 短桥 \([0.0015,0.004]\) 是与 Q50 重叠的独立交叉验证。

因此，当前 \(\kappa=4\) 基准的 \(Q_Z\) 全时间唯一零点问题已经关闭。该结论不自动推广到其它 Robin 参数或几何。下一阶段的核心任务是将此结果接入 \(h_T'<h_X'\) 的中间时间证明，并完成原来的 \(r(t)\) 单峰定理。

## 9. 中间时间窗口的 \(Q_{30}\) 单穿越区间证书

本节针对当前具体基准
\[
\kappa=4,\quad
w_1=25/66,\ D_1=25,\quad
w_2=4/33,\ D_2=1/4,\quad
D_R=1
\]
完成一个有限时间窗口认证。注意这里仍然是静态等阻力基准，但
\[
L_1\ne R_1,\qquad c=\lambda_2,\qquad L_1<c<R_1.
\]
因此必须使用本节前面修正后的 \(C_0\) 和留数公式，不能把 \(c\) 当成两个半谱的共同根。

取 \(N=30\)，定义
\[
Q_{30}(t)=
\sum_{2\le i<j\le30}
A_iA_j(\rho_j-\rho_i)^2e^{-(\rho_i+\rho_j)t}.
\]
用 70 位 mpmath 区间算术认证了 \(L_1,R_1,\lambda_2,R_2,\ldots,R_{30}\) 的根括区。右半根由严格递增的相位
\[
\theta_R(k)=k/2-\arctan(4/k),\qquad
\theta_R'(k)=1/2+\frac4{k^2+16}>0
\]
按相位单元编号。左半第一根由严格递增的左 Prüfer 相位认证。全区间第二根位于 \(L_1\) 与 \(R_1\) 之间，由中点 Dirichlet 分割的 Sturm–Liouville 交错确定其编号。

区间计算得到
\[
\begin{aligned}
k_{L_1}&\in[3.7726764010984087420047054886942,\\
&\qquad 3.7726764010984087420047056886943],\\
k_{R_1}&\in[4.5778594562068087296428282443474,\\
&\qquad 4.5778594562068087296428284443475],\\
\sqrt c&\in[3.9750626491471593943833679708026,\\
&\qquad 3.9750626491471593943833681708028].
\end{aligned}
\]
特别地
\[
c\in[15.80112306464483282514722870,\,
15.80112306464483282514723030].
\]
由归一化半谱乘积求得
\[
C_0\in[0.82836894199071981019641144,\,
0.82836894199071981019641583].
\]
留数上界常数可以取
\[
\boxed{B<265.424329828643.}
\]
所以
\[
|A_j|\le Bj,\qquad
\pi^2j^2\le\rho_j\le4\pi^2j^2
\]
对全部 \(j\ge2\) 成立。

### 9.1 有限和的严格单调性

对 \(Q_{30}\) 作区间求值，得到
\[
Q_{30}(0.010)\in
[-5525249.7078210,\,-5525249.7078209]<0,
\]
\[
Q_{30}(0.012)\in
[1955432.1495587,\,1955432.1495589]>0.
\]
将 \([0.010,0.012]\) 等分成 2000 个闭区间，并在每个区间上直接区间求值 \(Q_{30}'\)，取所有区间下端点的最小值，得到
\[
\boxed{
Q_{30}'(t)>1.0568131810\times10^9,
\qquad 0.010\le t\le0.012.
}
\tag{9.1}
\]
因此 \(Q_{30}\) 在整个窗口严格递增，且恰有一个零点。进一步，
\[
Q_{30}(0.01102)<0,\qquad
Q_{30}(0.01103)>0,
\]
故其唯一零点被严格括在
\[
\boxed{t_{30,*}\in(0.01102,0.01103).}
\tag{9.2}
\]

### 9.2 无限谱尾不能改变穿越次数

使用 (4.4)–(4.5)，在 \(\tau=0.010\) 处计算
\[
U_m(\tau),\quad U_{m,30}(\tau),\qquad m=0,1,2,3.
\]
每个 Gaussian 级数先求和至 \(j=160\)，再用严格递减的相邻项比值作几何尾界。得到
\[
\boxed{
\sup_{t\ge0.010}|Q_Z(t)-Q_{30}(t)|
<1.654\times10^{-25},
}
\tag{9.3}
\]
以及
\[
\boxed{
\sup_{t\ge0.010}|Q_Z'(t)-Q_{30}'(t)|
<6.354\times10^{-21}.
}
\tag{9.4}
\]
因此在整个窗口内
\[
\boxed{
Q_Z'(t)>1.0568131809\times10^9>0.
}
\tag{9.5}
\]
并且 \(Q_Z(0.010)<0<Q_Z(0.012)\)。更窄的括区也保留：
\[
Q_Z(0.01102)<0<Q_Z(0.01103).
\]
所以完整无限谱曲率分子在该窗口内恰有一个简单零点，且由负变正：
\[
\boxed{
\exists!\ z_*\in(0.01102,0.01103):
\quad Q_Z(z_*)=0,\quad Q_Z'(z_*)>0.
}
\tag{9.6}
\]
因为 \(f_Z>0\)，这也是 \((\log f_Z)''\) 在该窗口内唯一的零点，并且在此处由负变正。

本节单独给出 \([0.010,0.012]\) 的中间窗口证书。结合第 8 节列出的短时、负曲率外侧区间和晚时正性证书，当前 \(\kappa=4\) 静态打平基准的全时间唯一简单零点结论已经闭合。本结论尚未单独证明 \(h_T'<h_X'\) 或最终 \(r(t)\) 单峰性。

### 9.3 可复核文件

- robin_q30_interval_certificate.py：70 位区间算术，根括区、留数、\(Q_{30}\)、\(Q_{30}'\) 与 Gaussian 尾界。
- robin_q30_interval_certificate_input.json：精确几何、截断阶数、时间窗口和验收标准。
- robin_q30_interval_certificate_output.txt：原生运行输出。
- robin_q30_interval_certificate_provenance.json：源码、输入、输出 SHA-256 与运行环境。

运行环境为 Python 3.13.5、mpmath 1.3.0、Linux x86_64；没有重算全谱，认证器只求取 \(R_2,\ldots,R_{30}\) 和 \(L_1,R_1,\lambda_2\) 的窄区间。
