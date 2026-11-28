# Robin 反向风险率导数比较：精确一维加权积分归约

## 1. 目标

令
\[
T_c=X_c+Z_c,
\]
其中 \(X_c,Z_c\) 独立，CDF 和密度分别为 \(F_X,F_T\) 与 \(f_X,f_T\)。反向风险率为
\[
h_X(t)=\frac{f_X(t)}{F_X(t)},\qquad
h_T(t)=\frac{f_T(t)}{F_T(t)}.
\]
定义条件分布
\[
\pi_t(dz)=
\frac{F_X(t-z)f_Z(z)}{F_T(t)}
\mathbf 1_{(0,t)}(z)\,dz,
\]
并记
\[
q_t(z)=h_X(t-z),\qquad A_t(z)=-h_X'(t-z)>0.
\]
已有精确恒等式
\[
h_T'=\mathbb E_{\pi_t}[h_X'(t-Z)]
+\operatorname{Var}_{\pi_t}(h_X(t-Z)).
\]
目标是证明
\[
\boxed{h_X'(t)-h_T'(t)>0.}
\tag{1.1}
\]

本笔记给出一个新的等价表示。它不要求条件密度全局对数凹，也不要求 \(Z_c\) 的 score 曲率保持固定符号。

## 2. 密度导数比值恒等式

定义
\[
\psi_X(u)=\frac{f_X'(u)}{F_X(u)}.
\]
由于严格正光学距离给出 \(f_X(0+)=0\)，卷积密度满足
\[
f_T'(t)=(f_X'*f_Z)(t).
\]
另一方面，
\[
\begin{aligned}
\mathbb E_{\pi_t}[\psi_X(t-Z)]
&=\frac1{F_T(t)}
\int_0^t
F_X(t-z)f_Z(z)
\frac{f_X'(t-z)}{F_X(t-z)}\,dz\\
&=\frac{f_T'(t)}{F_T(t)}.
\end{aligned}
\]
故
\[
h_T'=
\mathbb E_{\pi_t}[\psi_X(t-Z)]-h_T(t)^2.
\]
又
\[
h_X'=\psi_X-h_X^2.
\]
相减得到
\[
\boxed{
h_X'(t)-h_T'(t)
=
h_T(t)^2-h_X(t)^2
+\psi_X(t)-\mathbb E_{\pi_t}[\psi_X(t-Z)].
}
\tag{2.1}
\]

这是条件方差恒等式的另一种精确形式。

### 2.1 一个简单但不全局适用的充分条件

如果 \(\psi_X\) 在 \((0,t)\) 单调不减，那么 \(Z_c>0\) 给出
\[
\mathbb E_{\pi_t}[\psi_X(t-Z)]
\leq\psi_X(t).
\]
而已有反向风险率序
\[
h_T(t)>h_X(t)>0
\]
于是推出
\[
h_X'(t)-h_T'(t)>0.
\]

因此
\[
\boxed{\psi_X'(u)\geq0\quad(u>0)}
\tag{2.2}
\]
是一个足以解决目标的充分条件。

但该条件不能全局成立。

## 3. 单调性障碍：由两个时间端点严格确定

记 \(d_X=\ell_L=7/22\)。已有高频 Robin 反演给出带可微余项的渐近：
\[
f_X(t)\sim C_Xt^{-3/2}
e^{-d_X^2/(4t)},
\]
\[
F_X(t)\sim \frac{4t^2}{d_X^2}f_X(t).
\]
因此
\[
\psi_X(t)=\frac{f_X'(t)}{F_X(t)}
\sim \frac{d_X^4}{16t^4},
\]
进而
\[
\boxed{
\psi_X'(t)\sim-\frac{d_X^4}{4t^5}<0,
\qquad t\downarrow0.
}
\tag{3.1}
\]

另一方面，设
\[
\eta=\lambda_3-c>0.
\]
最慢极点留数给出
\[
f_X(t)\sim B_Xe^{-\eta t},\qquad F_X(t)\longrightarrow1.
\]
所以
\[
\boxed{
\psi_X'(t)\sim
\eta^2B_Xe^{-\eta t}>0,
\qquad t\to\infty.
}
\tag{3.2}
\]

由 (3.1)–(3.2)，\(\psi_X'\) 至少发生一次符号变化。因此，直接证明 (2.2) 的全时间方案不成立。

这不是 \(h_T'<h_X'\) 的反例，而是排除了一个过强的充分条件。

## 4. 新的递增权重表示

令
\[
U_t=t-Z_c
\quad\text{在条件 }X_c+Z_c\leq t\text{ 下},
\]
并记它的 CDF 为 \(w_t\)。由条件密度直接得到
\[
\boxed{
w_t(u)=
\frac1{F_T(t)}
\int_0^u F_X(v)f_Z(t-v)\,dv,
\quad 0<u<t.
}
\tag{4.1}
\]
于是
\[
w_t(0)=0,\qquad w_t(t)=1,\qquad w_t'(u)>0.
\]

从 (2.1) 出发，
\[
\psi_X(t)-\psi_X(t-Z_c)
=\int_{t-Z_c}^t\psi_X'(u)\,du.
\]
交换期望与积分得到
\[
\psi_X(t)-\mathbb E_{\pi_t}\psi_X(t-Z_c)
=\int_0^t w_t(u)\psi_X'(u)\,du.
\tag{4.2}
\]

同样，
\[
h_T(t)-h_X(t)
=\mathbb E_{\pi_t}[h_X(t-Z_c)-h_X(t)]
=\int_0^t w_t(u)[-h_X'(u)]\,du.
\tag{4.3}
\]
将 (4.2)–(4.3) 代入 (2.1)，得到本轮核心恒等式：
\[
\boxed{
\begin{aligned}
h_X'(t)-h_T'(t)
=\int_0^t w_t(u)\Big[
&\psi_X'(u)\\
&+(h_T(t)+h_X(t))[-h_X'(u)]
\Big]\,du.
\end{aligned}
}
\tag{4.4}
\]

这里的权重 \(w_t\) 严格递增。因而负的短时 \(\psi_X'\) 并不是以均匀权重进入目标积分；它受到靠近原点时迅速消失的 CDF 权重抑制。

## 5. 等价的单积分表示

将 (4.4) 与 (4.1) 代入并交换积分次序，得到
\[
\boxed{
\begin{aligned}
h_X'(t)-h_T'(t)
=\frac1{F_T(t)}
\int_0^t F_X(v)f_Z(t-v)
\Big\{&
\psi_X(t)-\psi_X(v)\\
&+[h_T(t)+h_X(t)]
[h_X(v)-h_X(t)]
\Big\}\,dv.
\end{aligned}
}
\tag{5.1}
\]

定义
\[
\chi_t(v)=
\psi_X(v)-[h_T(t)+h_X(t)]h_X(v).
\]
则花括号正好是
\[
\chi_t(t)-\chi_t(v).
\]
所以
\[
\boxed{
h_X'-h_T'
=\mathbb E_{\omega_t}[\chi_t(t)-\chi_t(V)],
}
\tag{5.2}
\]
其中
\[
\omega_t(dv)=
\frac{F_X(v)f_Z(t-v)}{F_T(t)}
\mathbf1_{(0,t)}(v)\,dv
\]
是概率测度。

因此目标等价于
\[
\boxed{
\mathbb E_{\omega_t}[\chi_t(V)]<\chi_t(t).
}
\tag{5.3}
\]

这给出一个比直接控制 \(\operatorname{Var}_{\pi_t}(q_t)\) 更具体的研究对象：函数 \(\chi_t\) 与 Robin 卷积诱导测度 \(\omega_t\) 的比较。

## 6. 原点权重的解析抑制

固定 \(t>0\)。由 (4.1)，当 \(u\downarrow0\) 时，
\[
f_Z(t-u)=f_Z(t)+O(u).
\]
同时
\[
F_X(u)\asymp u^{1/2}e^{-d_X^2/(4u)}.
\]
端点 Laplace 估计给出
\[
\boxed{
w_t(u)=
O\left(u^{5/2}e^{-d_X^2/(4u)}\right),
\qquad u\downarrow0.
}
\tag{6.1}
\]
结合
\[
\psi_X'(u)=O(u^{-5}),
\]
可得
\[
w_t(u)\psi_X'(u)
=O\left(u^{-5/2}e^{-d_X^2/(4u)}\right),
\]
故 (4.4) 中的原点积分绝对收敛。

这解释了一个重要现象：虽然 \(\psi_X'\) 在原点附近具有强负奇异性，但它在目标积分中受到 \(F_X\) 所产生的光学指数抑制。

## 7. 当前可用的严格结论与下一步

本轮严格建立了：

1. \(h_X'-h_T'\) 的密度导数比值恒等式 (2.1)；
2. 全局单调性充分条件 (2.2)；
3. 该充分条件在短时必然失效、在晚时反向成立；
4. 递增条件权重表示 (4.4)；
5. 完全等价的单积分比较 (5.1)；
6. 原点奇异项受到显式光学指数抑制。

**尚未证明：**
\[
h_X'(t)-h_T'(t)>0
\]
在全部 \(0<t<t_*\) 上成立。

下一步不再试图证明 \(\psi_X\) 全局单调，而应研究
\[
\chi_t(t)-\mathbb E_{\omega_t}\chi_t(V)>0.
\]
重点是利用 \(F_X(v)\) 的小时间指数抑制，以及 \(f_Z(t-v)\) 的谱尾结构，控制 \(\chi_t\) 在早时间区域的负贡献。

## 8. 文献与计算状态

本轮复用了已有的 Robin 谱更新核、条件方差和首次通过时间文献审查；没有发现需要改变研究方向的新定理，也没有重新计算谱根。

本轮是解析推导，没有新的数值结果。所有先前已认证的谱曲率、根区间和谱尾证书继续沿用。
