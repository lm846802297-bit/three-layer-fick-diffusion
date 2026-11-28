# Robin 反向风险率导数比较：精确加权积分归约

## 1. 目标与记号

令 \(T_c=X_c+Z_c\)，其中 \(X_c,Z_c\) 独立。记 CDF 与密度为 \(F_X,F_T\) 和 \(f_X,f_T\)，反向风险率为
\[
h_X=\frac{f_X}{F_X},\qquad h_T=\frac{f_T}{F_T}.
\]
条件延迟分布为
\[
\pi_t(dz)=\frac{F_X(t-z)f_Z(z)}{F_T(t)}\,dz,\quad 0<z<t.
\]
已有
\[
h_T'=\mathbb E_{\pi_t}[h_X'(t-Z_c)]
+\operatorname{Var}_{\pi_t}(h_X(t-Z_c)).
\]
目标是证明 \(h_T'<h_X'\)。

## 2. 密度导数比值恒等式

定义
\[
\psi_X(u)=\frac{f_X'(u)}{F_X(u)}.
\]
严格正光学距离给出 \(f_X(0+)=0\)，故 \(f_T'=f_X'*f_Z\)。于是
\[
\mathbb E_{\pi_t}[\psi_X(t-Z_c)]
=\frac{f_T'(t)}{F_T(t)}
\]
并且
\[
h_T'=\mathbb E_{\pi_t}[\psi_X(t-Z_c)]-h_T^2,
\qquad h_X'=\psi_X-h_X^2.
\]
因此
\[
\boxed{
h_X'-h_T'
=h_T^2-h_X^2+\psi_X(t)-\mathbb E_{\pi_t}\psi_X(t-Z_c).
}
\tag{2.1}
\]

若 \(\psi_X\) 在 \((0,t)\) 单调不减，则 \(Z_c>0\) 和 \(h_T>h_X>0\) 立即推出目标不等式。这是一个充分条件。

## 3. 全局单调性充分条件不成立

令 \(d_X=\ell_L=7/22\)。短时高频反演给出
\[
f_X(t)\sim C_Xt^{-3/2}e^{-d_X^2/(4t)},\qquad
F_X(t)\sim\frac{4t^2}{d_X^2}f_X(t).
\]
所以
\[
\psi_X(t)\sim\frac{d_X^4}{16t^4},
\qquad
\boxed{\psi_X'(t)\sim-\frac{d_X^4}{4t^5}<0}
\quad(t\downarrow0).
\]
晚时首极点给出 \(f_X(t)\sim B_Xe^{-\eta t}\)、\(F_X(t)\to1\)，从而
\[
\boxed{\psi_X'(t)\sim\eta^2B_Xe^{-\eta t}>0}
\quad(t\to\infty).
\]
故不能用 \(\psi_X'\geq0\) 作为全时间证明。

## 4. 递增条件权重表示

令 \(U_t=t-Z_c\)，其分布取在条件 \(X_c+Z_c\leq t\) 下，并记 CDF 为 \(w_t\)。直接得到
\[
\boxed{
w_t(u)=\frac1{F_T(t)}
\int_0^u F_X(v)f_Z(t-v)\,dv,
\quad 0<u<t.
}
\tag{4.1}
\]
于是 \(w_t(0)=0,\ w_t(t)=1,\ w_t'(u)>0\)。

由
\[
\psi_X(t)-\psi_X(t-Z_c)=\int_{t-Z_c}^t\psi_X'(u)\,du
\]
以及
\[
h_T-h_X(t)=\int_0^t w_t(u)[-h_X'(u)]\,du,
\]
得到本轮核心恒等式
\[
\boxed{
h_X'(t)-h_T'(t)
=\int_0^t w_t(u)
\left\{
\psi_X'(u)+[h_T(t)+h_X(t)][-h_X'(u)]
\right\}du.
}
\tag{4.2}
\]
这是条件方差竞争的精确一维加权形式。

交换积分次序，还可写成
\[
\boxed{
\begin{aligned}
h_X'-h_T'
=\frac1{F_T(t)}\int_0^t F_X(v)f_Z(t-v)
\Big\{&
\psi_X(t)-\psi_X(v)\\
&+[h_T(t)+h_X(t)][h_X(v)-h_X(t)]
\Big\}\,dv.
\end{aligned}
}
\tag{4.3}
\]

## 5. 原点奇异项受到指数权重抑制

固定 \(t>0\)，由 \(f_Z(t-u)=f_Z(t)+O(u)\) 和
\[
F_X(u)\asymp u^{1/2}e^{-d_X^2/(4u)}
\]
得到
\[
\boxed{
w_t(u)=O\left(u^{5/2}e^{-d_X^2/(4u)}\right),
\quad u\downarrow0.
}
\]
而 \(\psi_X'(u)=O(u^{-5})\)，所以 \(w_t(u)\psi_X'(u)\) 绝对可积。短时负奇异性虽使简单单调性判据失效，但在真实条件权重下被光学指数抑制。

## 6. 研究状态

本轮严格建立了 (2.1)、(4.2)、(4.3)，并排除了 \(\psi_X\) 全局单调这一过强充分条件。尚未证明 (4.2) 的积分严格为正。

下一步应控制
\[
\psi_X'(u)+[h_T(t)+h_X(t)][-h_X'(u)]
\]
的负区间，并利用 \(w_t(u)\) 在原点的指数小量估计其总贡献。已有 \(f_Z\) 全时间唯一曲率零点定理可用于进一步控制 \(f_Z(t-v)\) 的变化；不需要重新计算谱根。


## 7. 利用已证明的 \(f_Z\) 曲率定理控制延迟权重

令
\[
\ell_Z(s)=(\log f_Z)'(s).
\]
对固定 \(\kappa=4\) 基准，已有全时间定理给出
\[
\ell_Z'(s)>0,\qquad s>z_*,
\]
而晚时首谱项给出
\[
\ell_Z(s)\longrightarrow-\rho_2.
\]
故任意 \(\varepsilon>z_*\) 都满足
\[
\ell_Z(\varepsilon)\leq\ell_Z(s)<-\rho_2,\qquad s\geq\varepsilon.
\]
令 \(M_\varepsilon=-\ell_Z(\varepsilon)\)。对 \(t\geq2\varepsilon\) 和 \(0\leq v\leq\varepsilon\)，积分 score 得
\[
\boxed{
e^{\rho_2v}\leq\frac{f_Z(t-v)}{f_Z(t)}
\leq e^{M_\varepsilon v}.
}
\tag{7.1}
\]

### 7.1 一个有限标量验收判据

假设存在 \(u_0<\varepsilon\)，使
\[
\psi_X'(u)<0\quad(0<u<u_0),\qquad
\psi_X'(u)\geq0\quad(u\geq u_0).
\tag{7.2}
\]
记
\[
\begin{aligned}
\mathcal D(\varepsilon)
={}&[\psi_X(2\varepsilon)-\psi_X(\varepsilon)]
\int_0^\varepsilon F_X(v)e^{\rho_2v}\,dv\\
&-\int_0^{u_0}
F_X(v)e^{M_\varepsilon v}
[\psi_X(v)-\psi_X(u_0)]\,dv.
\end{aligned}
\tag{7.3}
\]

由 (4.2)，令
\[
B_t(u)=\psi_X'(u)+[h_T(t)+h_X(t)][-h_X'(u)].
\]
因为第二项严格为正，(7.2) 蕴含 \(B_t(u)\geq0\) 对 \(u\geq u_0\) 成立，而在 \(u<u_0\) 上 \(B_t(u)\geq\psi_X'(u)\)。

对 \(t\geq2\varepsilon\)，保留正区间 \([\varepsilon,2\varepsilon]\)，利用 \(w_t\) 递增：
\[
\int_\varepsilon^{2\varepsilon}w_t(u)B_t(u)\,du
\geq w_t(\varepsilon)[\psi_X(2\varepsilon)-\psi_X(\varepsilon)].
\]
负区间利用 (4.1) 交换积分次序，并使用 (7.1) 上界；正区间权重 \(w_t(\varepsilon)\) 使用 (7.1) 下界。于是
\[
\boxed{
h_X'(t)-h_T'(t)
\geq \frac{f_Z(t)}{F_T(t)}\mathcal D(\varepsilon),
\qquad t\geq2\varepsilon.
}
\tag{7.4}
\]
所以只要 (7.2) 和
\[
\boxed{\mathcal D(\varepsilon)>0}
\tag{7.5}
\]
得到严格认证，就关闭全部 \(t\geq2\varepsilon\) 的反向风险率导数比较。

### 7.2 下一步的固定验收任务

当前工作候选为
\[
\varepsilon=0.035,\qquad 2\varepsilon=0.070.
\]
该值来自有限谱探索，用于选择后续认证窗口；它不是已经通过区间算术验收的常数，也不构成定理前提。

下一步只做三项：
1. 对 \(\psi_X'\) 建立短时负性与 \(u\geq u_0\) 非负性的严格分区证书，并严格括定 \(u_0\)；
2. 以 \(\varepsilon=0.035\) 为首个候选，使用已有根区间、有限谱和解析尾界认证 \(\mathcal D(\varepsilon)>0\)；若裕量不足，只在同一标量判据内调整 \(\varepsilon\)；
3. 独立关闭 \(0<t\leq0.070\) 的剩余区间。

这样，原来的全时间 \(h_T'<h_X'\) 问题就被拆成一个短时间问题、一个 \(\psi_X'\) 符号问题和一个有限标量积分判据。