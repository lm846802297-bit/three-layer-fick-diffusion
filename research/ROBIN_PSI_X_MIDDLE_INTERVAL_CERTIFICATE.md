# Robin \(X_c\) 中间时间区间的严格负性证书

## 1. 定理

固定静态等阻力三层 Robin 基准
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区宽度 \(1/2\)、扩散率 \(1\)，Robin 参数 \(\kappa=4\)。令
\[
\eta_i=\lambda_{i+1}-c,\qquad c=\lambda_2,
\]
并写
\[
f_X(t)=\sum_{i\ge2}B_i e^{-\eta_i t},\qquad
F_X(t)=1-\sum_{i\ge2}\frac{B_i}{\eta_i}e^{-\eta_i t}.
\]
定义
\[
N_X(t)=f_X''(t)F_X(t)-f_X'(t)f_X(t),\qquad
\psi_X(t)=\frac{f_X'(t)}{F_X(t)}.
\]
本证书证明
\[
\boxed{N_X(t)<-0.01122<0,\qquad
0.0015\le t\le0.010.}
\tag{1.1}
\]
由于 \(F_X(t)>0\)，这等价于
\[
\boxed{\psi_X'(t)<0,\qquad
0.0015\le t\le0.010.}
\]

## 2. 精确留数与根区间

全区间特征函数 \(\Delta_G(k)\) 由三层传递矩阵及 Robin 端点条件定义，右半特征函数为
\[
Y_R(k)=\cos\frac{k}{2}+\frac4k\sin\frac{k}{2}.
\]
归一化特征函数为 \(P_G(z)=\Delta_G(\sqrt{-z})/24\)、\(P_R(z)=Y_R(\sqrt{-z})/3\)。对 \(n=i+1\)，精确留数公式为
\[
B_i=C
\frac{P_R(-\lambda_{i+1})}{1-\lambda_{i+1}/R_1}
\frac{(1-\lambda_{i+1}/\lambda_1)(1-\lambda_{i+1}/c)}
{P_G'(-\lambda_{i+1})},
\]
其中
\[
C=\frac{1-c/R_1}{P_R(-c)}
\frac{cP_G'(-c)}{1-c/\lambda_1},\qquad
P_G'(-\lambda_n)=-\frac{\Delta_G'(k_n)}{48k_n}.
\]
本轮对前 82 个全区间波数构造区间根括，半宽为 \(10^{-30}\)，并用区间 Prüfer 相位计数认证相应 Sturm 节点数恰为 \(0,1,\ldots,81\)。因此所用根确为全谱前 82 根，没有漏根或错序。

首个右半根 \(R_1\) 直接复用此前 Q30 根区间证书，不重新计算。由区间特征函数及其导数得到前 79 个留数 \(B_i\) 和速率 \(\eta_i\) 的包含区间。

## 3. 无限谱尾的解析界

### 3.1 速率下界

全算子二次型为
\[
\mathfrak q[u]=\int_0^1D(x)|u'(x)|^2dx
+\kappa(|u(0)|^2+|u(1)|^2).
\]
在两个界面释放连续性约束并去掉正 Robin 边界项，得到三个独立 Neumann 区间算子。由 min-max 原理，原算子的特征值不小于该 Neumann 直和的特征值。其光学总长度为
\[
L_{\rm opt}=\frac{w_1}{\sqrt{D_1}}+
\frac{w_2}{\sqrt{D_2}}+\frac{1/2}{1}
=\frac9{11}.
\]
Neumann 直和计数函数满足
\[
N_{\rm Neu}(\Lambda)\le
\frac{L_{\rm opt}}{\pi}\sqrt{\Lambda}+3,
\]
故
\[
\lambda_n\ge
\frac{\pi^2(n-3)^2}{L_{\rm opt}^2}.
\]
取 \(n=i+1\)，利用
\[
\frac{121\pi^2}{81}>14.74,
\]
以及在 \(i=81\) 的正裕量和随后递增的差分，得到
\[
\boxed{\eta_i=\lambda_{i+1}-c\ge14i^2,
\qquad i\ge81.}
\tag{3.1}
\]

### 3.2 留数增长界

由全区间 Sturm–Liouville 的 Lagrange 恒等式，
\[
|P_G'(-\lambda_n)|=
\frac{\int_0^1y_n(x)^2dx}{24|y_n(1)|}.
\]
当 \(k_n\ge10\) 时，第一层显式解给出
\[
\int_0^1y_n^2dx>0.029.
\]
对缩放状态 \((y,p/k)\)，三层传递矩阵无穷范数依次不超过 \(6,3,2\)，故 \(|y_n(1)|\le36\)。因此
\[
|P_G'(-\lambda_n)|>\frac{0.029}{24\cdot36}.
\]
同时 \( |P_R(-\lambda_n)|\le1.4/3\)、\(|C|<1\)，且对 \(i\ge81\)，\(\lambda_{i+1}\ge2\max(R_1,\lambda_1,c)\)。代入留数公式可得
\[
|B_i|<3\times10^4\lambda_{i+1}.
\]
由 Dirichlet 试探子空间，
\[
\lambda_{i+1}\le25\pi^2(i+1)^2
\le100\pi^2i^2,
\]
于是
\[
\boxed{|B_i|\le3\times10^7i^2,
\qquad i\ge81.}
\tag{3.2}
\]

### 3.3 Gaussian 型尾项

对 \(m=0,1,2\) 及 \(t\ge\tau=0.0015\)，
\[
\left|\sum_{i\ge81}B_i(-\eta_i)^m e^{-\eta_i t}\right|
\le3\times10^7(100\pi^2)^m
\sum_{i\ge81}i^{2m+2}e^{-14i^2\tau}.
\]
利用首项与递减几何比，得到
\[
\begin{aligned}
T_0&<2.961\times10^{-49},\\
T_1&<1.919\times10^{-42},\\
T_2&<1.244\times10^{-35}.
\end{aligned}
\]
对 \(F_X\) 的尾项另有
\[
\left|\sum_{i\ge81}\frac{B_i}{\eta_i}e^{-\eta_i t}\right|
<3.221\times10^{-54}.
\]
Taylor 余项所需的 \(m=0,\ldots,11\) 统一绝对导数界均在原生输出中记录。

## 4. 时间区间认证

使用前 79 个留数 \(i=2,\ldots,80\)，并把上述无限谱尾作为对称区间误差加入。对
\[
N_X=f_X''F_X-f_X'f_X
\]
在每个时间盒中心构造八阶 Taylor 多项式。其导数系数由
\[
\begin{aligned}
N_X^{(n)}={}&\sum_{k=0}^n\binom nk
f_X^{(k+2)}F_X^{(n-k)}\\
&-\sum_{k=0}^n\binom nk
f_X^{(k+1)}f_X^{(n-k)}
\end{aligned}
\]
精确计算；九阶余项用全区间绝对导数界控制。

验收参数：

- 区间运算精度 80 位；根与留数的中间计算精度 90 位；
- Taylor 阶数 8；
- 时间盒 426 个；
- 最大盒宽小于 \(2.0\times10^{-5}\)；
- 全局九阶导数余项界小于 \(1.5931\times10^{43}\)；
- 每盒 Taylor 余项小于 \(4.2982\times10^{-8}\)。

所有 426 个盒均通过严格负性检查，汇总上界为
\[
\sup_{[0.0015,0.010]}N_X(t)
\le-0.0112285080871\ldots<-0.01122.
\]
故 (1.1) 得证。

## 5. 与已完成区间的连接

此前的晚时间区间证书已经证明
\[
\begin{aligned}
N_X(t)&<0 &&(0.010\le t\le0.0144),\\
N_X'(t)&>0 &&(0.0144\le t\le0.0145),\\
N_X(t)&>0 &&(t\ge0.0145).
\end{aligned}
\]
因此本证书将严格负性向前推进到 \(t=0.0015\)，从而严格得到
\[
\boxed{\psi_X'(t)<0\quad(0.0015\le t\le0.0144),}
\]
且在 \((0.0144,0.0145)\) 内恰好一次由负变正，随后保持正值。

**本证书没有覆盖 \(0<t<0.0015\)。** 先前对话中曾提出一个基于反射路径的 \(N_X\) 短时估计，但尚无与之对应的可复核源码、输入和原生输出，因此不将其计入本定理。当前可严格确认的结论是 \(\psi_X'\) 在 \([0.0015,\infty)\) 上恰好一次变号。

## 6. 可复核文件

- robin_psi_x_middle_interval_certificate.py：区间根、留数、谱尾及 Taylor-box 认证源码。
- robin_psi_x_middle_interval_certificate_input.json：完整输入与验收参数。
- robin_psi_x_middle_interval_certificate_output.txt：本机原生输出。
- robin_psi_x_middle_interval_certificate_provenance.json：源码、输入、输出哈希与运行环境。

运行环境：Python 3.13.5、mpmath 1.3.0、NumPy 2.3.5、SciPy 1.17.0，Linux x86_64，5 个可见 CPU 核心；未使用 GPU。
