# 三层 Fick 扩散：三个反超可以一起移到任意晚

## 定理（局部结论）

仍用 \(\Delta=A_L-A_R\)、\(\Delta(t)=d-\sum_{n\ge1}c_n e^{-\lambda_nt}\)，其中 \(d=\Delta_\infty\)。在单位区间、中点源、理想界面及等 Robin 壁面条件下，存在另一个临界点：左壁快层 \(D=4\)，厚度 \(h\simeq0.20448868701134\)；左内层 \(D=B\simeq0.6583335027765\)，厚度 \(0.5-h\)；右半层 \(D=1\)；\(\kappa\simeq7.6171345835013\)。这里 \(d=c_1=c_2=0\)，\(c_3<0\)。

临界点附近有开邻域，所有参数均在 \(0<t\le0.3\) 左壁领先，之后至多有三次反超。对任意给定的观测截止时间，邻域中存在参数，使**恰好三次**反超全部发生在该截止之后。此为原晚期双反超临界点之外的第二个结果，不是全族上限。

## 临界参数的认证

在 \(d=0\) 截面，\(B=(0.5-h)/(0.5-h/4)\)。从两壁分别向中点传播状态 \((y,q)=(1,\kappa)\)，每层厚度 \(w\)、阻力 \(r=w/D\)，用矩阵
\[
\begin{pmatrix}
\cos X & r\,\sin X/X\\
-\lambda w\,\sin X/X & \cos X
\end{pmatrix},\quad X=\sqrt{\lambda wr}.
\]
中点两侧为 \((a,b),(c,e)\)，全区间本征条件为 \(W=ae+bc=0\)。用 Krawczyk 夹区和区间雅可比认证以下四个方程：
\[
W(\lambda_1)=0,\qquad a(\lambda_1)-c(\lambda_1)=0,\qquad a(\lambda_2)=c(\lambda_2)=0.
\]
第一模两壁等振幅；第二模在中点有它唯一的节点。谱索引通过分层 Sturm 节点数核对；对第二模使用每个半区间的第一 Robin–Dirichlet 零点，避免把界面节点重复计算。第三模系数 \(c_3\simeq-0.0123902982552\)。

### 局部坐标

\(d\) 对 \(B\) 的导数非零。就打平面而言，先用 \(W(\lambda_1)=0\) 和 \(c(\lambda_2)=0\) 消去两个简单根；四元雅可比非退化意味着 \(a(\lambda_1)-c(\lambda_1)\) 和 \(a(\lambda_R)\) 在 \((h,\kappa)\) 中独立，\(\lambda_R\) 为右半区间 Robin–Dirichlet 本征值。前一差值与 \(c_1\) 的微分成非零比例。第二模处，对 \(W=0\) 微分，令 \(Q=a(\lambda_R)\)，则全区间本征根处的中点值的微分为
\[
{\rm d}a_{\rm eig}=\frac{b\,c_\lambda}{W_\lambda}\,{\rm d}Q .
\]
比例非零，因为两侧节点处通量非零，两类本征根都简单；且第二模在两个壁面异号，故 \(c_2\) 的微分与 \({\rm d}a_{\rm eig}\) 成非零比例。于是 \((d,c_1,c_2)\) 是三维局部坐标。

## 晚期下限：构造恰好三个指定的渐近交叉时刻

设 \(T\) 足够大，用局部坐标取 \(d=\alpha_0 e^{-\lambda_3 T}\) 和 \(c_j=\alpha_j e^{-(\lambda_3-\lambda_j)T}\)（\(j=1,2\)），右端本征值取临界值。对有界的 \(\tau\)，
\[
e^{\lambda_3 T}\Delta(T+\tau)
\longrightarrow \alpha_0-\alpha_1 e^{-\lambda_1\tau}
-\alpha_2 e^{-\lambda_2\tau}-c_3 e^{-\lambda_3\tau}.
\]
参数偏差指数小，本征值移位乘以 \(T\) 仍趋零，高阶谱尾也一致趋零。用前三项插值，使极限函数在任意三个预先确定的不同 \(\tau\) 处为零；不同指数构成 Chebyshev 系统（反复除以最慢指数并用 Rolle 定理），保证这三个零点简单且两侧异号，故三次过零继承到足够大的有限 \(T\)。

## 早期及晚期上限

短时反射余项严格给出 \(0<t\le0.001\) 的正号；用 140 个有编号的模态、区间 Taylor 包络与统一谱尾界，核验 \(0.000999\le t\le0.3\) 始终严格正。证书在临界点上运行，严格裕量和参数连续性给出开邻域；未给显式最大邻域。

令 \(G=e^{\lambda_3t}\Delta\)。依次对 \(G\) 施加
\[
(\partial_t-(\lambda_3-\lambda_2))(\partial_t-(\lambda_3-\lambda_1))(\partial_t-\lambda_3).
\]
三次加权微分恰好消去 \(d,c_1,c_2\)，只留下第三模与更高模。用区间认证对 \(t\ge0.3\) 有
\[
\sum_{n\ge4}|c_n|\lambda_n(\lambda_n-\lambda_1)(\lambda_n-\lambda_2)
 e^{-(\lambda_n-\lambda_3)t}
<|c_3|\lambda_3(\lambda_3-\lambda_1)(\lambda_3-\lambda_2).
\]
只需在 \(t=0.3\) 认证，因为左端递减；高阶和的误差界在参数邻域一致有效。结果是三次加权微分后的函数没有零点。每次加权微分可写成乘正指数、求导、再乘正指数，应用 Rolle 定理三遍给出晚期至多三个零点（含相切也不超过三个不同零点）。与早期正号及上述下限构造合并，得到恰好三个任意晚的反超。

运行：`OPENBLAS_NUM_THREADS=1 python late_triple.py`（Python、NumPy、SciPy、mpmath）；`fold.py` 仅用于定位根。短时估计推导沿用附带 EXACTLY_THREE.md。
