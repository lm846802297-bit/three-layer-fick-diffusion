# 全时间唯一零点定理：静态打平 Robin 基准 \(\kappa=4\)

## 定理

考虑静态三层 Robin 扩散基准
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区长度 \(1/2\)、扩散率 \(D_R=1\)，两端 Robin 参数
\[
\kappa=4.
\]
令 \(L_j,R_j\) 分别为左右半区 Robin–Dirichlet 谱，令
\[
c=\lambda_2,\qquad L_1<c<R_1,\qquad
\rho_j=R_j-c,\quad j\geq2.
\]
定义归一化延迟密度
\[
f_Z(t)=\sum_{j\geq2}A_je^{-\rho_jt}>0
\]
及其曲率分子
\[
Q_Z(t)=f_Z(t)f_Z''(t)-f_Z'(t)^2
=\sum_{2\leq i<j}A_iA_j(\rho_j-\rho_i)^2
e^{-(\rho_i+\rho_j)t}.
\]
则存在唯一
\[
\boxed{z_*\in(0.01102,0.01103)}
\]
使得
\[
\boxed{
Q_Z(t)
\begin{cases}
<0,&0<t<z_*,\\
=0,&t=z_*,\\
>0,&t>z_*.
\end{cases}}
\]
等价地，\((\log f_Z)''\) 在正半轴恰有一个简单零点，并由负变正。

## 证明

全时间轴由四个互相衔接的认证区间覆盖。

### I. 短时间：\(0<t\leq0.0015\)

使用精确双曲传递矩阵分解
\[
\Phi_Z(s)=K e^{-\delta q}H(q)\frac{N_L(q)}{N_R(q)},
\qquad q=\sqrt{s-c},\quad \delta=\frac2{11},
\]
其中
\[
K=\frac{11C_0L_1}{2R_1}.
\]
因 \(L_1\ne R_1\)，此处必须保留首半谱根剥离因子 \(L_1/R_1\)。

将主项写为三个正卷积比的组合，并将界面反射项按最短额外光学路径
\[
d_0=\frac13
\]
分离。反射相对误差具有指数因子
\[
e^{-\alpha/t},\qquad
\alpha=\frac{d_0^2-\delta^2}{4}
=0.0195133149678604\ldots.
\]
显式 Bromwich 估计控制反射误差及其前两阶时间导数。70 位区间常数认证给出
\[
(\log f_0)''<-\frac{0.0139}{t^3},
\]
\[
\left|(\log(1+\varepsilon))''\right|<48410.72,
\qquad f_Z=f_0(1+\varepsilon).
\]
因此
\[
\boxed{
(\log f_Z)''(t)<-4.0701\times10^6<0,
\quad 0<t\leq0.0015.
}
\]
完整推导见短时反射路径证书。

### II. 负曲率区间：\(0.0015\leq t\leq0.010\)

取 \(N=50\)，构造有限谱和
\[
Q_{50}(t)=\sum_{2\leq i<j\leq50}
A_iA_j(\rho_j-\rho_i)^2e^{-(\rho_i+\rho_j)t}.
\]
使用 70 位区间根括区和 300 个 Taylor-box 时间子区间，得到
\[
\sup Q_{50}
\leq-5\,517\,850.2501393.
\]
由 \(|A_j|\leq Bj\)、\(\pi^2j^2<\rho_j<4\pi^2j^2\) 导出的 Gaussian 双谱尾界为
\[
\sup_{t\geq0.0015}|Q_Z-Q_{50}|<61.579377296.
\]
故
\[
\boxed{
Q_Z(t)<-5\,517\,788.6707<0,
\quad 0.0015\leq t\leq0.010.
}
\]

### III. 唯一穿越窗口：\(0.010\leq t\leq0.012\)

已有 \(N=30\) 的 70 位区间证书证明
\[
Q_Z'(t)>1.0568131809\times10^9>0
\]
在整个窗口成立，并且
\[
Q_Z(0.01102)<0<Q_Z(0.01103).
\]
因此窗口内恰有一个简单零点，且
\[
z_*\in(0.01102,0.01103).
\]

### IV. 晚时间：\(t\geq0.012\)

首对谱项满足
\[
C_{23}=A_2A_3(\rho_3-\rho_2)^2>0,
\qquad \lambda_{23}=\rho_2+\rho_3
\]
且其指数率严格小于所有其它成对谱率。定义
\[
E(t)=
\sum_{(i,j)\ne(2,3)}
\frac{|A_iA_j|(\rho_j-\rho_i)^2}{C_{23}}
e^{-[(\rho_i+\rho_j)-\lambda_{23}]t}.
\]
在 \(t=0.012\) 的有限谱与无限谱尾区间估计给出
\[
E(0.012)<0.762157<1.
\]
由于 \(E(t)\) 严格递减，
\[
Q_Z(t)\geq C_{23}e^{-\lambda_{23}t}[1-E(t)]>0,
\qquad t\geq0.012.
\]

四段区间覆盖整个 \(t>0\)，且唯一穿越窗口中的零点严格单调穿越，定理得证。 \(\square\)

## 可复核证书

| 区间 | 证书 | 核心结果 |
|---|---|---|
| \((0,0.0015]\) | ROBIN_SHORT_TIME_REFLECTION_CERTIFICATE.md | 反射路径解析余项，严格负性 |
| \([0.0015,0.010]\) | ROBIN_CURVATURE_SIGN_EXTENSION.md | Q50 Taylor-box + Gaussian 尾界 |
| \([0.010,0.012]\) | robin_q30_interval_certificate.py | Q30 的 \(C^1\) 单穿越证书 |
| \([0.012,\infty)\) | robin_q50_taylor_box_300.py | 首对谱项支配，\(E(0.012)<1\) |

额外的 Q30 短桥 \([0.0015,0.004]\) 与 Q50 区间重叠，提供独立交叉验证。

每份证书的精确输入、原生输出、运行环境和 SHA-256 均记录在各自的 Provenance 文件中。本定理没有重新计算全谱；它组合已有根区间、解析谱尾界和分区间认证结果。

## 结论范围

本定理严格针对上述固定三层几何及 \(\kappa=4\) 基准。它不自动推广到其它 Robin 参数或几何，也不单独推出
\[
h_T'(t)<h_X'(t)
\]
或最终的 \(r(t)\) 单峰性。后两者仍是下一阶段的独立证明目标。
