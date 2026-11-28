# \(\psi_X'\) 的晚时间单次穿越区间证书：\(\kappa=4\) 静态打平基准

## 1. 认证范围

令
\[
\psi_X(t)=\frac{f_X'(t)}{F_X(t)},\qquad
N_X(t)=f_X''(t)F_X(t)-f_X'(t)f_X(t).
\]
由于 \(F_X(t)>0\)，
\[
\operatorname{sgn}\psi_X'(t)=\operatorname{sgn}N_X(t).
\]
本证书只覆盖 \(t\ge0.010\)，不把尚未处理的 \((0,0.010)\) 纳入结论。

<badge color="success">区间认证通过</badge>

\[
\boxed{
\begin{aligned}
N_X(t)&< -456.20,
&&0.010\le t\le0.0144,\\
N_X'(t)&>1.5569\times10^7,
&&0.0144\le t\le0.0145,\\
N_X(t)&>1056.56,
&&0.0145\le t\le0.020.
\end{aligned}}
\tag{1.1}
\]
因此 \(N_X\) 在 \((0.0144,0.0145)\) 恰有一个简单零点。

## 2. 谱根与留数的区间构造

沿用前一份精确留数推导：
\[
B_i=C
\frac{P_R(-\lambda_{i+1})}{1-\lambda_{i+1}/R_1}
\frac{(1-\lambda_{i+1}/\lambda_1)(1-\lambda_{i+1}/c)}
{P_G'(-\lambda_{i+1})}.
\]
本轮使用 70 位区间算术包围 102 个全区间特征根和 101 个右半区间特征根；全区间根由 Sturm 节点数 \(0,1,\ldots,101\) 逐个索引，右半根位于
\[
(2j-1)\pi<k_{R,j}<2j\pi.
\]
根区间半宽为 \(10^{-35}\)。由特征函数和导数的区间传递矩阵，计算前 100 个 \(B_i,\eta_i\) 的包含区间。

## 3. 无限谱尾的解析界

对全区间特征函数使用 Lagrange 恒等式：
\[
|P_G'(-\lambda_n)|
=\frac{\int_0^1 y_n(x)^2\,dx}{24|y_n(1)|}.
\]
在 \(k_n\ge10\) 时，第一层显式解给出
\[
\int_0^1y_n^2dx
\ge \frac12\left(\frac{w_1}{2}-\frac{5}{4k_n}\right)
-\frac{16w_1}{25k_n^2}
>0.029.
\]
对缩放状态 \((y,p/k)\)，三层传递矩阵的无穷范数乘积不超过
\[
6\times3\times2=36,
\]
故 \(|y_n(1)|\le36\)，从而
\[
|P_G'(-\lambda_n)|>\frac{0.029}{24\cdot36}.
\]
结合
\[
|P_R(-\lambda_n)|\le\frac{1.4}{3},\qquad
|C|<1,
\]
以及 \(\lambda_n\ge2\max(R_1,\lambda_1,c)\)，留数公式推出
\[
|B_i|\le 25211\,\lambda_{i+1}.
\]
由 Dirichlet 上界和 Neumann 下界，
\[
\lambda_{i+1}\le100\pi^2i^2,\qquad
\eta_i\ge2.4i^2,\quad i\ge101.
\]
因此采用稍宽松的统一界
\[
\boxed{|B_i|\le3\times10^7i^2,\quad
\eta_i\le100\pi^2i^2,\quad
\eta_i\ge2.4i^2\quad(i\ge101).}
\tag{3.1}
\]
对 \(m=0,1,2,3\)，尾项满足
\[
\left|\sum_{i>100}B_i(-\eta_i)^m e^{-\eta_it}\right|
\le
3\times10^7(100\pi^2)^m
\sum_{i\ge101}i^{2m+2}e^{-2.4i^2t}.
\tag{3.2}
\]
右端用首项与递减几何比控制；在 \(t\ge0.010\) 时远小于上述区间符号裕量。

## 4. 晚时间首两模态支配

在 \(t_0=0.020\)，区间计算得到
\[
\begin{aligned}
-\bigl(f_X'\bigr)_{\{2,3\}}(t_0)
&\in[1356.2070172520368,\ 1356.2070172520385],\\
\sum_{i\ge4}|B_i|\eta_i e^{-\eta_it_0}
&<292.874,\\[2mm]
\bigl(f_X''\bigr)_{\{2,3\}}(t_0)
&>100310.60,\\
\sum_{i\ge4}|B_i|\eta_i^2 e^{-\eta_it_0}
&<66462.63.
\end{aligned}
\]
因为 \(\eta_i>\eta_3\) 对所有 \(i\ge4\)，高阶绝对尾与首两模态正贡献之比随 \(t\) 严格下降。因此
\[
\boxed{f_X'(t)<0,\qquad f_X''(t)>0,\quad t\ge0.020.}
\tag{4.1}
\]
结合 \(f_X,F_X>0\)，有
\[
N_X=f_X''F_X-f_X'f_X>0,\qquad t\ge0.020.
\]

## 5. 本轮定理

综合 (1.1) 与 (4.1)，得到
\[
\boxed{
\psi_X'(t)<0\quad(0.010\le t\le0.0144),
}
\]
\[
\boxed{
\psi_X'(t)>0\quad(t\ge0.0145),
}
\]
并且 \(\psi_X'\) 在 \((0.0144,0.0145)\) 内恰好一次由负变正。

因此，**\(\psi_X\) 在 \([0.010,\infty)\) 上恰有一个严格极小点**，其位置满足
\[
u_0\in(0.0144,0.0145).
\]

尚未关闭的区间只有
\[
\boxed{0<t<0.010.}
\]
短时渐近已给出
\[
\psi_X'(t)\sim-\frac{(7/22)^4}{4t^5}<0,
\]
但还需要一个显式可微余项界，将该渐近与 \(t=0.010\) 接合。故目前不能把唯一极小点结论推广到整个正半轴。
