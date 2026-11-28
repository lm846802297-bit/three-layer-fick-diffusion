# Robin \(X_c\) 密度导数比值 \(\psi_X=f_X'/F_X\)：精确留数与单极小点探索

## 1. 本轮命题

沿用固定三层静态打平基准
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区宽度 \(1/2\)、扩散率 \(1\)，Robin 参数 \(\kappa=4\)。令
\[
c=\lambda_2,\qquad
\eta_i=\lambda_{i+1}-c,\qquad
\rho_i=R_i-c,\quad i\ge2.
\]
\(X_c\) 的 Laplace 变换为
\[
\widehat f_X(s)=
\prod_{i\ge2}\frac{1+s/\rho_i}{1+s/\eta_i}.
\]
定义
\[
\psi_X(t)=\frac{f_X'(t)}{F_X(t)},\qquad
N_X(t)=f_X''(t)F_X(t)-f_X'(t)f_X(t).
\]
则精确地
\[
\boxed{\psi_X'(t)=\frac{N_X(t)}{F_X(t)^2}.}
\tag{1.1}
\]
因此 \(\psi_X'\) 的符号完全由 \(N_X\) 决定。并且
\[
\boxed{N_X'(t)=f_X'''(t)F_X(t)-f_X'(t)^2.}
\tag{1.2}
\]
式 (1.2) 是认证唯一穿越窗口时可直接使用的导数公式。

## 2. 通过特征函数得到精确留数

定义全区间特征函数 \(\Delta_G(k)\)：从左端初值
\[
(y,p)=(1,\kappa)
\]
出发，依次乘以两层左半区和均匀右半区的传递矩阵，记末端状态为 \((y_1,p_1)\)，并令
\[
\Delta_G(k)=p_1+\kappa y_1.
\]
右半区特征函数为
\[
Y_R(k)=\cos\frac{k}{2}+\frac{\kappa}{k}\sin\frac{k}{2}.
\]
归一化特征乘积写成
\[
P_G(z)=\frac{\Delta_G(\sqrt{-z})}{24}
=\prod_{n\ge1}\left(1+\frac{z}{\lambda_n}\right),
\]
\[
P_R(z)=\frac{Y_R(\sqrt{-z})}{3}
=\prod_{j\ge1}\left(1+\frac{z}{R_j}\right).
\]
这里 \(\Delta_G(0)=24\)、\(Y_R(0)=3\)。

剥离首个右半根及全谱前两个根：
\[
P_R^{(1)}(z)=\frac{P_R(z)}{1+z/R_1},
\qquad
P_G^{(12)}(z)=
\frac{P_G(z)}
{(1+z/\lambda_1)(1+z/c)}.
\]
于是有精确恒等式
\[
\boxed{
\widehat f_X(s)=
\frac{P_R^{(1)}(s-c)}{P_G^{(12)}(s-c)}
\frac{P_G^{(12)}(-c)}{P_R^{(1)}(-c)}.
}
\tag{2.1}
\]

令 \(n=i+1\)。在极点 \(s=-\eta_i\) 处取留数，得到
\[
B_i=
\frac{P_G^{(12)}(-c)}{P_R^{(1)}(-c)}
\frac{P_R^{(1)}(-\lambda_n)}
{(P_G^{(12)})'(-\lambda_n)}.
\]
消去剥离因子后，
\[
\boxed{
B_i=C
\frac{P_R(-\lambda_{i+1})}{1-\lambda_{i+1}/R_1}
\frac{(1-\lambda_{i+1}/\lambda_1)
      (1-\lambda_{i+1}/c)}
     {P_G'(-\lambda_{i+1})},
}
\tag{2.2}
\]
其中
\[
\boxed{
C=
\frac{1-c/R_1}{P_R(-c)}
\frac{cP_G'(-c)}{1-c/\lambda_1}.
}
\tag{2.3}
\]
对 \(k_n=\sqrt{\lambda_n}\)，导数可直接从传递矩阵计算：
\[
P_G'(-\lambda_n)
=-\frac{\Delta_G'(k_n)}{48k_n}.
\tag{2.4}
\]
这避免了用有限谱乘积近似无限乘积留数所产生的系统性偏差。

在正时间上，
\[
f_X(t)=\sum_{i\ge2}B_i e^{-\eta_i t},
\qquad
F_X(t)=1-\sum_{i\ge2}\frac{B_i}{\eta_i}e^{-\eta_i t}.
\tag{2.5}
\]
由全区间 Sturm–Liouville 的 Lagrange 恒等式，在全谱根 \(\lambda_n=k_n^2\) 处有
\[
|P_G'(-\lambda_n)|
=\frac{\int_0^1 y_n(x)^2\,dx}{24|y_n(1)|}.
\]
对缩放状态 \((y,p/k)\)，每层传递矩阵在高频下具有一致有界的矩阵范数，故 \(|y_n(1)|\le C\)。而第一层内
\[
y_n(x)=\cos(k_nx/5)+\frac{\kappa}{5k_n}\sin(k_nx/5),
\]
所以 \(\int_0^1y_n^2dx\ge c_0>0\) 对充分大的 \(n\) 成立。于是 \(|P_G'(-\lambda_n)|\ge c_1>0\)。代入 (2.2)，并使用 \(\lambda_n=O(n^2)\)，得到多项式留数界
\[
|B_i|\le C_2 i^2.
\]
因此 (2.5) 的留数级数及其任意固定阶时间导数在每个 \(t\ge\tau>0\) 上正常收敛。

## 3. 两端符号已由解析渐近确定

高频传递矩阵给出 \(X_c\) 的光学长度
\[
d_X=\ell_L=\frac7{22}.
\]
因此
\[
f_X(t)\sim C_Xt^{-3/2}e^{-d_X^2/(4t)},
\qquad
F_X(t)\sim\frac{4t^2}{d_X^2}f_X(t).
\]
由此
\[
\boxed{
\psi_X'(t)\sim-\frac{d_X^4}{4t^5}<0,
\qquad t\downarrow0.
}
\tag{3.1}
\]
晚时首极点支配给出
\[
f_X(t)\sim B_2e^{-\eta_2t},\qquad B_2>0,
\]
从而
\[
\boxed{
\psi_X'(t)\sim\eta_2^2B_2e^{-\eta_2t}>0,
\qquad t\to\infty.
}
\tag{3.2}
\]
所以至少存在一个正时间零点。唯一性仍需证明。

## 4. 精确留数的高精度探索

使用 60 位 mpmath，先由传递矩阵隔离全谱前 62 个波数和右半谱前 61 个波数，再用 (2.2) 直接计算精确特征函数留数。探索结果为
\[
\begin{aligned}
\lambda_1&=5.85153416383962693\ldots,\\
c=\lambda_2&=15.80112306464483283\ldots,\\
\lambda_3&=53.22857301669635033\ldots,\\
R_1&=20.95679720078209853\ldots,\\
R_2&=103.5096693904747410\ldots,\\
B_2&=41.95278570484199704\ldots,\\
B_3&=55.18793423410773466\ldots.
\end{aligned}
\]
求解 \(N_X(t)=0\) 得到候选值
\[
\boxed{
u_0^{\rm num}=0.0144307180068361996801201055\ldots
}
\tag{4.1}
\]
并且
\[
N_X(0.0144)<0,\qquad N_X(0.0145)>0.
\]
在候选根处的探索值为
\[
N_X'(u_0^{\rm num})
\approx1.5886947042\times10^7>0.
\]
因此数值证据支持一个简单的负转正穿越。

根位置的截断稳定性也很强：10 模给出 \(0.0144307180798\ldots\)，15 模给出 \(0.0144307180068\ldots\)，20 模及以上在当前精度下给出相同的根位置。

### 4.1 一个重要的高频警告

在 \(t=0.0015\)，10、15、20 模截断给出的 \(\psi_X'\) 仍然为正；到 30 模才转为负，40–60 模稳定在
\[
\psi_X'(0.0015)\approx-2.9317701905\times10^{11}.
\]
这说明小时间符号由高阶模的集体抵消决定。不能用低模的符号或低阶截断来证明 (3.1) 与中间区间的连接。

本节全部数值均为高精度探索值，不是区间包含证书，也不把 (4.1) 宣称为已认证根区间。

## 5. 文献边界

Z. Szabó, “On Monotonic Behaviour of Relative Increments of Unimodal Distributions,” *Acta Mathematica Academiae Paedagogicae Nyíregyháziensis* **22** (2006), 113–119, [全文](http://real.mtak.hu/186874/1/amapn22_12.pdf)，研究单峰密度的相对增量函数和普通风险率的单调阶段，并给出基于 \(f,f',f''\) 的判据。将变量反射 \(X\mapsto-X\) 可把 \(f_X/F_X\) 视作反射变量的普通风险率；但该文没有直接给出
\[
\left(\frac{f_X'}{F_X}\right)'
\]
恰好一次变号的定理。

Keilson 关于首次通过时间密度的对数凹/对数凸结果同样针对单个 passage-time density，不直接控制这里的密度导数比值。当前尚无可直接复用的定理证明 \(\psi_X\) 的唯一极小点。

## 6. 当前证明状态与唯一剩余区间

后续已新增 [\(\psi_X'\) 晚时间单次穿越区间证书](ROBIN_PSI_X_LATE_SINGLE_CROSSING_CERTIFICATE.md)。它以 70 位区间算术和解析谱尾界证明：
\[
\begin{aligned}
N_X&<0 &&(0.010\le t\le0.0144),\\
N_X'&>0 &&(0.0144\le t\le0.0145),\\
N_X&>0 &&(t\ge0.0145).
\end{aligned}
\]
因此现在已经严格证明
\[
\boxed{\psi_X'\text{ 在 }[0.010,\infty)\text{ 上恰好一次由负变正}.}
\]
零点被括在 \((0.0144,0.0145)\)。

**唯一剩余的符号区间是**
\[
\boxed{0<t<0.010.}
\]
短时高频渐近给出
\[
\psi_X'(t)\sim-\frac{(7/22)^4}{4t^5}<0,
\]
但尚需将可微短时反射余项显式控制到 \(t=0.010\)，或用一段有限谱区间证书与该渐近接合。

完成这一短时接合后，才能把 \(\psi_X'\) 的全时间单次变号性质代入前一轮的 \(\mathcal D(\varepsilon)\) 加权尾项判据。当前不提前宣称 \(\psi_X\) 的全时间唯一极小点定理已经完成。