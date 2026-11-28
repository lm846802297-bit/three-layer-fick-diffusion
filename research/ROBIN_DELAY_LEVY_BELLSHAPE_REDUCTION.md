# Robin 延迟谱乘积的 Lévy 表示与 bell-shaped 判据：静态打平基准的反例

## 1. 对象

沿用指数倾斜后的延迟变量
\[
\widehat\nu_c(s)=
\prod_{i\ge2}\frac{1+s/\sigma_i}{1+s/\rho_i},
\qquad
\rho_i=R_i-c,\quad \sigma_i=L_i-c,
\]
其中
\[
0<\rho_i\le \sigma_i,\qquad
\rho_i,\sigma_i\asymp i^2.
\]
严格光学长度差给出
\[
\prod_{i\ge2}\frac{\rho_i}{\sigma_i}=0,
\]
故 \(Z_c\) 无零点原子。

## 2. 精确 Lévy–Khintchine 表示

逐项使用
\[
\log\frac{1+s/\rho_i}{1+s/\sigma_i}
=
\int_{\rho_i}^{\sigma_i}\frac{s}{u(u+s)}\,du
=
\int_0^\infty(1-e^{-sx})
\left(\int_{\rho_i}^{\sigma_i}e^{-ux}\,du\right)dx,
\]
得到
\[
\boxed{
-\log\widehat\nu_c(s)
=
\int_0^\infty(1-e^{-sx})\Pi_c(x)\,dx,
}
\]
其中
\[
\boxed{
\Pi_c(x)=
\sum_{i\ge2}\int_{\rho_i}^{\sigma_i}e^{-ux}\,du
=
\sum_{i\ge2}\frac{e^{-\rho_i x}-e^{-\sigma_i x}}{x}.
}
\tag{2.1}
\]
因此 \(Z_c\) 是无漂移次级过程在单位时间的分布。

对每个 \(n\ge0\)，
\[
(-1)^n\Pi_c^{(n)}(x)
=
\sum_{i\ge2}\int_{\rho_i}^{\sigma_i}u^ne^{-ux}\,du>0.
\]
所以
\[
\boxed{\Pi_c\text{ 是严格完全单调的 Lévy 密度。}}
\]
这将 \(Z_c\) 归入具有完全单调 Lévy 密度的 Bondesson 类。

定义谱区间重叠计数
\[
\boxed{
m(u)=\sum_{i\ge2}\mathbf1_{[\rho_i,\sigma_i)}(u).
}
\]
则
\[
\Pi_c(x)=\int_0^\infty e^{-ux}m(u)\,du.
\]
由 \(\rho_i\asymp i^2\)，
\[
\int_1^\infty\frac{m(u)}{u^3}\,du
\le \frac12\sum_{i\ge2}\rho_i^{-2}<\infty,
\qquad
\int_1^\infty\frac{m(u)}{u^2}\,du
\le\sum_{i\ge2}\rho_i^{-1}<\infty.
\tag{2.2}
\]

## 3. Kwaśnicki–Simon 表示与必要条件

取特征函数
\[
\Phi(\xi)=\mathbb E e^{-i\xi Z_c}=\widehat\nu_c(i\xi).
\]
由 Lévy 表示，
\[
\log\Phi(\xi)
=
\int_0^\infty
\left(\frac1{u+i\xi}-\frac1u\right)m(u)\,du.
\]
将 \(u>1\) 的一阶项分离：
\[
\begin{aligned}
\log\Phi(\xi)
={}&-ib\xi\\
&+\int_0^\infty
\left[
\frac1{u+i\xi}
-\left(\frac1u-\frac{i\xi}{u^2}\right)\mathbf1_{\{u>1\}}
\right]m(u)\,du,
\end{aligned}
\]
其中
\[
b=\int_1^\infty\frac{m(u)}{u^2}\,du<\infty.
\]
所以 Kwaśnicki–Simon 表示中的 level-crossing 函数正是
\[
\varphi(u)=
\begin{cases}
m(u),&u>0,\\
0,&u<0.
\end{cases}
\]

正则条件成立，因为 \(|\Phi(\xi)|\le1\) 且 \(\Phi(0)=1\)；可积条件由 (2.2) 给出。正光学间隙还给出 Fourier 变换的高频衰减，故 \(Z_c\) 具有光滑密度。

Kwaśnicki–Simon 的完整刻画不只是一个充分判据：每个 weakly bell-shaped 函数的 Fourier 表示都必须满足整数水平穿越条件，而且表示函数由 Fourier 变换唯一恢复。参见：

- M. Kwaśnicki, *A new class of bell-shaped functions*, Transactions of the AMS 373 (2020), 349–382, [DOI](https://doi.org/10.1090/tran/7825).
- M. Kwaśnicki and T. Simon, *Characterisation of the class of bell-shaped functions*, Mathematische Zeitschrift 301 (2022), 3207–3230, [DOI](https://doi.org/10.1007/s00209-022-02997-7), 尤其是 Theorem 1.3 及其前面的 Theorem 1.1。

因此，如果 \(m(u)-k\) 对某个整数 \(k\) 出现“负—正—负”符号序列，则 \(f_Z\) 不可能是 weakly bell-shaped，也就不可能是严格 bell-shaped 密度。

## 4. 整数水平穿越的精确谱条件

由于 \((\rho_i)\)、\((\sigma_i)\) 均严格递增，在固定 \(u\) 处活跃的谱区间指标构成连续整数段。于是
\[
m(u)\ge k+1
\iff
\exists i\ge2:\quad \rho_{i+k}\le u<\sigma_i.
\tag{4.1}
\]
令
\[
d_{i,k}:=\sqrt{L_i}-\sqrt{R_{i+k}}.
\]
则 Kwaśnicki–Simon 的 level-crossing 条件等价于
\[
\boxed{
d_{i,k}>0,\ i<j
\quad\Longrightarrow\quad
d_{j,k}\ge0,
\qquad k\ge1.
}
\tag{4.2}
\]
注意，等号只对应 \(m(u)=k\) 的零平台，不构成负向穿越。真正禁止的是 \(d_{i,k}\) 先为正、后来又为负。

## 5. 现有静态打平 Robin 基准：\(k=12\) 的严格反例

取已经归档的静态打平基准：
\[
w_1=\frac{25}{66},\quad D_1=25,
\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区
\[
D_R=1,\qquad \kappa=4.
\]
其光学长度为
\[
\ell_L=\frac7{22},\qquad \ell_R=\frac12.
\]
本节只使用已经存在的模型参数，不重新求解整个谱。

### 5.1 四个半谱根的区间括定

使用单调 Prüfer 相位方程，得到以下波数根区间：
\[
\begin{aligned}
\sqrt{L_{21}}&\in[205.1617247434,\ 205.1617247438],\\
\sqrt{L_{22}}&\in[209.5748600934,\ 209.5748600939],\\
\sqrt{R_{33}}&\in[204.2426865656,\ 204.2426865661],\\
\sqrt{R_{34}}&\in[210.5247035095,\ 210.5247035100].
\end{aligned}
\tag{5.1}
\]
相位残差在各括区两端严格异号；左半相位导数严格为正，右半相位导数也严格为正，因此这些括区分别对应所标记的第 21、22、33、34 个根。

此外，对同一组括区直接代入两层传递矩阵特征函数，端点仍严格异号。这是与相位残差不同表达式的独立核对。

区间运算给出
\[
\boxed{
\sqrt{L_{21}}-\sqrt{R_{33}}
\in[0.9190381773,\ 0.9190381782]>0,
}
\]
以及
\[
\boxed{
\sqrt{L_{22}}-\sqrt{R_{34}}
\in[-0.9498434166,\ -0.9498434155]<0.
}
\tag{5.2}
\]
因此
\[
\boxed{
\rho_{33}<\sigma_{21}<\sigma_{22}<\rho_{34}.
}
\tag{5.3}
\]
这里 \(\sigma_i=L_i-c\)、\(\rho_i=R_i-c\)，公共倾斜参数 \(c\) 在比较中完全抵消。

### 5.2 直接违反第 12 个整数水平条件

由 (5.3)，可选取
\[
u_+\in(\rho_{33},\sigma_{21}),
\qquad
u_-\in(\sigma_{22},\rho_{34}),
\]
且
\[
u_+<\sigma_{21}<\sigma_{22}<u_-.
\]

在 \(u_+\) 处，对所有 \(i=21,\ldots,33\)，均有
\[
\rho_i\le\rho_{33}<u_+<\sigma_{21}\le\sigma_i.
\]
因此至少有 13 个谱区间同时活跃：
\[
m(u_+)\ge13,
\qquad
m(u_+)-12>0.
\]

在 \(u_-\) 处，所有 \(i\le22\) 的区间已经退出，而所有 \(i\ge34\) 的区间尚未进入。因此至多只有 \(i=23,\ldots,33\) 这 11 个区间可能活跃：
\[
m(u_-)\le11,
\qquad
m(u_-)-12<0.
\]

同时，在 \(u_0<\rho_2\) 时 \(m(u_0)=0\)，故
\[
m(u_0)-12<0.
\]
于是严格得到
\[
\boxed{
m(u_0)-12<0,\qquad
m(u_+)-12>0,\qquad
m(u_-)-12<0,
\quad u_0<u_+<u_-.
}
\tag{5.4}
\]

这已经是第 12 个整数水平的两次符号变化，不满足 Kwaśnicki–Simon 条件。

### 5.3 结论

\[
\boxed{
\text{该静态打平 Robin 基准的 }f_Z
\text{ 不是 bell-shaped 密度。}
}
\tag{5.5}
\]
这是一个解析区间证书，而不是浮点扫描结论。

## 6. 对先前高频有限化论证的修正

此前曾尝试使用
\[
\sqrt{L_{i+1}}-\sqrt{L_i}\longrightarrow\frac{\pi}{\ell_L}
\]
来推出所有充分高阶相邻谱间距的统一下界。对于扩散系数在内部界面发生跳变的分层介质，这一步没有依据：Weyl 渐近控制的是平均谱计数，不自动控制每一个相邻波数间距。

在当前基准中，忽略 \(O(k^{-1})\) 的 Robin 边界修正，左半特征函数的高频主项为
\[
\begin{aligned}
Y_L(k)
&=\cos\left(\frac{5k}{66}\right)
  \cos\left(\frac{16k}{66}\right)
 -10\sin\left(\frac{5k}{66}\right)
  \sin\left(\frac{16k}{66}\right)
 +O(k^{-1})\\
&=\frac12\left[
11\cos\left(\frac{21k}{66}\right)
-9\cos\left(\frac{11k}{66}\right)
\right]+O(k^{-1}).
\end{aligned}
\tag{6.1}
\]
它是两个相位频率的组合，而不是单一余弦。故不能仅由总光学长度推出相邻左半波数间距逐项收敛到 \(\pi/\ell_L\)。

因此，原先基于统一相邻谱间距的 \(I_0,K_0\) 有限化证明不适用于当前分层基准。更重要的是，(5.4) 已经直接否定该基准的 bell-shape 条件，不必再试图通过高阶谱尾来补救。

## 7. 对原始曲率问题的影响

本结果只排除一种证明工具：
\[
f_Z\text{ 全局 bell-shaped}.
\]
它不意味着
\[
(\log f_Z)''
\]
没有可利用的局部符号结构，也不意味着最终的 \(r(t)\) 单峰性不成立。

因此，后续不再以 bell-shaped 密度为前提，也不再从 \(f_Z''\) 的两个拐点推导中央对数凹区间。应回到已经建立的精确恒等式
\[
h_T'-h_X'
=
A_t(0)+
\int_0^tW_t(z)(\log f_Z(z))''\,dz,
\]
直接研究谱留数和条件分布对这个加权积分的控制。

## 8. 可复核归档

新增证书文件：

- robin_bellshape_level12_obstruction.py：60 位 mpmath 区间算术；相位括根及独立传递矩阵端点符号核对。
- robin_bellshape_level12_obstruction_input.json：精确几何、目标整数水平及四个根括区。
- robin_bellshape_level12_obstruction_output.txt：原生区间运行结果。
- robin_bellshape_level12_obstruction_provenance.json：源码、输入、输出 SHA-256 与运行环境。

本轮没有重新计算已有扩散谱的全体根；只对四个已有高阶根位置进行了窄区间认证。当前静态打平基准的 bell-shape 检验已经关闭，后续研究回到加权 score 曲率问题。
