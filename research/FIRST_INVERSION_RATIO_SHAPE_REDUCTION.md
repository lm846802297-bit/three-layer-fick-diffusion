# 首模逆序：导数比值的精确卷积归约

## 1. 目标量

令
\[
F_L(t)=\frac{A_L(t)}{A_L(\infty)},\qquad
F_R(t)=\frac{A_R(t)}{A_R(\infty)},\qquad
\gamma=\frac{A_L(\infty)}{A_R(\infty)}>1,
\]
并定义
\[
G_L=F_L+\frac{F_L'}{\lambda_1},\qquad
G_R=F_R+\frac{F_R'}{\lambda_1}.
\]
则
\[
J=(\partial_t+\lambda_1)\Delta
=\lambda_1A_L(\infty)\left(G_L-\frac{G_R}{\gamma}\right).
\]
因此，令
\[
Q(t)=\frac{G_R(t)}{G_L(t)},
\]
有精确符号等价
\[
J(t)<0\quad\Longleftrightarrow\quad Q(t)>\gamma.
\]

目标是证明 \(Q\) 至多有一个最大值。若成立，集合 \(\{t>0:Q(t)>\gamma\}\) 至多是一个区间，从而 \(J\) 至多有两个正时间零点。

## 2. 半谱乘积的精确分解

在首模唯一逆序条件下写
\[
a=L_1,\qquad c=\lambda_2,\qquad b=R_1,
\qquad a<c<b,
\]
并沿用高阶半谱排序
\[
L_i\ge R_i\ge\lambda_{i+1},\qquad i\ge2.
\]
由捕获概率的精确 Laplace 变换及 \(G_L(0)=G_R(0)=0\)，得到
\[
\widehat{G_L'}(s)
=\prod_{i\ge1}\frac{1+s/R_i}{1+s/\lambda_{i+1}},
\]
\[
\widehat{G_R'}(s)
=\prod_{i\ge1}\frac{1+s/L_i}{1+s/\lambda_{i+1}}.
\]

定义正常高阶半谱核
\[
K(s)=\prod_{i\ge2}\frac{1+s/R_i}{1+s/\lambda_{i+1}},
\]
以及正等待分布 \(\nu\)：
\[
\widehat\nu(s)=
\prod_{i\ge2}\frac{1+s/L_i}{1+s/R_i}.
\]
每个因子均为概率变换：
\[
\frac{1+s/R_i}{1+s/\lambda_{i+1}}
=
\frac{\lambda_{i+1}}{R_i}
+\left(1-\frac{\lambda_{i+1}}{R_i}\right)
\frac{\lambda_{i+1}}{s+\lambda_{i+1}},
\]
\[
\frac{1+s/L_i}{1+s/R_i}
=
\frac{R_i}{L_i}
+\left(1-\frac{R_i}{L_i}\right)
\frac{R_i}{s+R_i}.
\]
故 \(K\) 是某个非负等待变量 \(X\) 的 Laplace 变换，\(\nu\) 是独立非负延迟的概率分布。

令 \(f_X\) 为 \(X\) 的密度，并定义
\[
k(t)=(e^{-c\,\cdot}*f_X)(t)
=\int_0^t e^{-c(t-u)}f_X(u)\,du.
\]
这正是静态打平时间域证明中使用的公共核；已有结果证明 \(k(t)>0\)、\(k(0+)=0\)，且 \(k\) 在 \(t>0\) 严格对数凹。

## 3. 新的精确时间域恒等式

首模因子分别满足
\[
\frac{1+s/b}{1+s/c}
=\frac cb\left(1+\frac{b-c}{s+c}\right),
\qquad
\frac{1+s/a}{1+s/c}
=\frac ca\left(1-\frac{c-a}{s+c}\right).
\]
利用 \(k'=f_X-ck\)，以及卷积与微分可交换，得到
\[
\boxed{
G_L'(t)=\frac cb\,[k'(t)+bk(t)] ,
}
\]
\[
\boxed{
G_R'(t)=\frac ca\,
\big[(\nu*k)'(t)+a(\nu*k)(t)\big]
=\frac ca\,[\nu*(k'+ak)](t).
}
\]

特别地，
\[
G_L'(t)=\frac cb
\left[f_X(t)+(b-c)k(t)\right]>0.
\]
积分并使用 \(k(0)=0\) 还得到
\[
G_L(t)=\frac cb
\left[k(t)+b\int_0^t k(u)\,du\right]>0,
\]
\[
G_R(t)=\frac ca
\left[(\nu*k)(t)+a\int_0^t(\nu*k)(u)\,du\right]>0.
\]
因此 \(Q(t)>0\)，且 \(G_L\) 可作为严格正的权重函数。

最关键的导数比值现在具有如下形式：
\[
\boxed{
r(t):=\frac{G_R'(t)}{G_L'(t)}
=\frac ba\,
\frac{\nu*(k'+ak)}{k'+bk}.
}
\]
这把原来的两个完整捕获响应之比，化成了一个正常半谱正延迟 \(\nu\)、一个已证明严格对数凹的公共核 \(k\)，以及两个一阶算子 \(D+a,D+b\) 的组合。

## 4. 单峰传递引理

由于 \(G_L(0)=G_R(0)=0\)，
\[
Q(t)=
\frac{\int_0^t r(u)G_L'(u)\,du}
{\int_0^t G_L'(u)\,du},
\]
所以 \(Q(t)\) 是 \(r\) 的历史加权平均。直接微分得
\[
Q'(t)=\frac{G_L'(t)}{G_L(t)}[r(t)-Q(t)].
\]
再令 \(D(t)=r(t)-Q(t)\)，可得
\[
\boxed{
\frac{d}{dt}\{G_L(t)D(t)\}=G_L(t)r'(t).
}
\]

由此得到一个可直接用于本模型的充分条件：

**单峰传递引理。** 若 \(r\) 在其正值区间上严格对数凹，且正值区间是一个从零开始的区间；在该区间之后 \(r\le0\)，则 \(Q\) 至多有一个最大值。

证明：对数凹性保证 \(r\) 在正值区间至多先增后减。若 \(r\) 递增，则它不小于此前的加权平均 \(Q\)，故 \(Q'\ge0\)。进入 \(r\) 的递减段后，\(G_L(r-Q)\) 单调不增，因此 \(r-Q\) 至多由正变负一次。此后若 \(r\le0\)，而 \(Q>0\)，则 \(Q'<0\)。所以 \(Q\) 至多有一个最大值。证毕。

因此，当前全域证明的核心已精确缩小为证明
\[
r(t)=\frac ba\frac{\nu*(k'+ak)}{k'+bk}
\]
在其正值区间上严格对数凹，并且正值区间只有一个连通分支。

## 5. 指数倾斜后的 CDF 表示：正值区间已经封闭

令
\[
M_X=K(-c),\qquad M_Z=\widehat\nu(-c),
\]
并定义指数倾斜概率测度
\[
dP_{X_c}(x)=\frac{e^{cx}}{M_X}\,dP_X(x),\qquad
dP_{Z_c}(z)=\frac{e^{cz}}{M_Z}\,d\nu(z).
\]
取 \(X_c,Z_c\) 独立，记 \(T_c=X_c+Z_c\)，其 CDF 分别为 \(F_X,F_T\)，密度为 \(f_X^c,f_T^c\)。为避免与未倾斜密度混淆，以下写
\[
h_X=\frac{f_X^c}{F_X},\qquad h_T=\frac{f_T^c}{F_T}.
\]
它们是反向风险率（reversed hazard rates）。

由指数倾斜的定义，
\[
\boxed{
k(t)=M_Xe^{-ct}F_X(t),\qquad
(\nu*k)(t)=M_XM_Ze^{-ct}F_T(t).
}
\]
因此
\[
k'+ak=M_Xe^{-ct}F_X(t)\,[h_X(t)-(c-a)],
\]
\[
(\nu*k)'+a(\nu*k)
=M_XM_Ze^{-ct}F_T(t)\,[h_T(t)-(c-a)],
\]
而
\[
k'+bk=M_Xe^{-ct}F_X(t)\,[h_X(t)+(b-c)]>0.
\]

### 5.1 倾斜后的阶段分布

对每个 \(i\ge2\)，\(X_c\) 的第 \(i\) 个独立阶段仍是零膨胀指数分布：
\[
\begin{array}{c|c|c}
&\text{零原子概率}&\text{非零时的指数速率}\\ \hline
X_c&
\dfrac{\lambda_{i+1}-c}{R_i-c}&
\lambda_{i+1}-c\\[6pt]
Z_c&
\dfrac{R_i-c}{L_i-c}&
R_i-c
\end{array}
\]
对应的非零概率分别为
\[
\frac{R_i-\lambda_{i+1}}{R_i-c},
\qquad
\frac{L_i-R_i}{L_i-c}.
\]
所有速率严格为正，因为
\[
c=\lambda_2<\lambda_3\le R_i\le L_i.
\]

任一零膨胀指数变量都可写成
\[
Y_i=(E_i-\tau_i)_+,\qquad E_i\sim\operatorname{Exp}(r_i),
\]
其中 \(\tau_i=-r_i^{-1}\log q_i\)，\(q_i\) 是该阶段的非零概率。有限阶段和的 CDF 是
\[
P\!\left(\sum_i(E_i-\tau_i)_+\le t\right).
\]
由于 \(\sum_i(E_i-\tau_i)_+\) 关于 \(E_i\) 是凸函数，其次水平集在 \((E,t)\) 空间中为凸集；指数联合密度是对数凹函数。由 Prékopa 定理，有限阶段和的 CDF 对数凹。又由 Weyl 渐近，两个无限阶段和的期望均收敛，有限部分和单调趋于 \(X_c,T_c\)，故取极限仍有
\[
\boxed{
F_X(t)\ \text{和}\ F_T(t)\ \text{均对数凹}.
}
\]

严格光学长度为正保证 \(X_c\) 没有零原子：其各阶段零原子概率趋于
\[
\left(\frac{\ell_R}{\ell_L+\ell_R}\right)^2<1,
\]
故无限乘积为零。于是 \(F_X(0)=F_T(0)=0\)，而 \(F_X,F_T>0\) 对所有 \(t>0\) 成立。热核/谱乘积的高频衰减给出正时间上的连续密度。因而
\[
h_X'(t)\le0,\qquad h_T'(t)\le0.
\]
此外
\[
h_T(0+)=+\infty,\qquad h_T(\infty)=0.
\]
第一式由 \(\log F_T(t)\to-\infty\) 和其凹性推出；第二式由 \(F_T(t)\to1\) 及 \(h_T\) 非增推出。

### 5.2 \(r\) 的正值集合是一个初始区间

令
\[
\alpha=c-a>0,\qquad \beta=b-c>0.
\]
代入第 3 节的精确卷积恒等式，得到
\[
\boxed{
r(t)=\frac ba\,M_Z\,
\frac{F_T(t)}{F_X(t)}
\frac{h_T(t)-\alpha}{h_X(t)+\beta}.
}
\]
分母严格为正，所以
\[
r(t)>0\quad\Longleftrightarrow\quad h_T(t)>\alpha.
\]
由于 \(h_T\) 连续、非增，且从 \(+\infty\) 降至 0，存在有限 \(t_*>0\)，使得
\[
\boxed{
r(t)>0\quad(0<t<t_*),\qquad
r(t)\le0\quad(t\ge t_*).
}
\]
若 \(h_T\) 在阈值处没有平台，则 \(t_*\) 是唯一解 \(h_T(t_*)=\alpha\)；即使存在平台，正值集合仍然只有一个连通分支。

因此，上一轮单峰传递引理的“正值区间从零开始且只有一个连通分支”这一前提已经证明，不再是开放条件。

### 5.3 一个额外的反向风险率排序

由
\[
\frac{(\nu*k)(t)}{k(t)}
=\int_{[0,t]}\frac{k(t-z)}{k(t)}\,\nu(dz)
\]
及 \(k\) 的对数凹性，对每个固定 \(z\ge0\)，
\[
\frac{d}{dt}\log\frac{k(t-z)}{k(t)}
=(\log k)'(t-z)-(\log k)'(t)\ge0.
\]
积分域也随 \(t\) 增大，故 \((\nu*k)/k\) 非减。结合倾斜表示，
\[
\boxed{
\frac{F_T(t)}{F_X(t)}\ \text{非减},
\qquad h_T(t)\ge h_X(t).
}
\]
这是 \(X_c\le_{\rm rh}T_c\) 的反向风险率序。它是后续曲率估计可直接使用的额外约束。

### 6.1 条件分布下的精确矩恒等式

定义
\[
\pi_t(dz)=\frac{F_X(t-z)\,\mathbf 1_{[0,t]}(z)\,\nu_c(dz)}{F_T(t)},
\]
其中 \(\nu_c=P_{Z_c}\)。这正是 \(Z_c\) 在条件 \(X_c+Z_c\le t\) 下的条件分布。令
\[
q_t(z)=h_X(t-z).
\]
由 \(F_T(t)=\int F_X(t-z)\nu_c(dz)\) 得
\[
h_T(t)=\mathbb E_{\pi_t}[q_t].
\]
对参数 \(t\) 求导时，权重的对数导数为 \(q_t-h_T\)。因此
\[
\boxed{
h_T'=\mathbb E_{\pi_t}[h_X'(t-Z)]
+\operatorname{Var}_{\pi_t}(h_X(t-Z)).
}
\]
这与 \(F_T\) 的对数凹性 \(h_T'\le0\) 一起给出严格可用的方差约束
\[
\boxed{
\operatorname{Var}_{\pi_t}(h_X(t-Z))
\le-\mathbb E_{\pi_t}[h_X'(t-Z)].
}
\]
再求一次导数，得到
\[
\boxed{
h_T''=
\mathbb E_{\pi_t}[h_X''(t-Z)]
+3\operatorname{Cov}_{\pi_t}
 \!\left(h_X(t-Z),h_X'(t-Z)\right)
+\kappa_{3,\pi_t}\!\left(h_X(t-Z)\right),
}
\]
其中 \(\kappa_{3,\pi_t}\) 是三阶中心矩。这些恒等式在正光学距离下成立；\(F_X(0+)=f_X^c(0+)=0\)，故移动积分端点不产生边界项。

特别地，
\[
h_T'-h_X'
=
\mathbb E_{\pi_t}[h_X'(t-Z)-h_X'(t)]
+\operatorname{Var}_{\pi_t}(h_X(t-Z)).
\]
所以，虽然已经证明 \(h_T\ge h_X\)，要进一步证明 \(F_T/F_X\) 对数凹，仍须控制上式中“平均曲率差”与条件方差之间的竞争。单凭 CDF 对数凹给出的方差上界还不足以确定这个差的符号。

### 5.4 两个倾斜分布的同速率 Bernoulli--指数表示

记
\[
\alpha_i=\lambda_{i+1}-c,\qquad
\rho_i=R_i-c,\qquad
\sigma_i=L_i-c.
\]
在首模唯一逆序及高阶排序条件下，
\[
0<\alpha_i\le \rho_i\le \sigma_i.
\]
第 5.1 节的倾斜阶段变换可以改写为
\[
\begin{aligned}
\mathbb E e^{-sX_{c,i}}
&=\frac{1+s/\rho_i}{1+s/\alpha_i}
 =\frac{\alpha_i}{\rho_i}
  +\left(1-\frac{\alpha_i}{\rho_i}\right)
    \frac{\alpha_i}{s+\alpha_i},\\
\mathbb E e^{-sZ_{c,i}}
&=\frac{1+s/\sigma_i}{1+s/\rho_i}
 =\frac{\rho_i}{\sigma_i}
  +\left(1-\frac{\rho_i}{\sigma_i}\right)
    \frac{\rho_i}{s+\rho_i}.
\end{aligned}
\]
因此可以在同一概率空间上取独立 \(E_i\sim\operatorname{Exp}(\alpha_i)\)，以及 Bernoulli 指示量 \(B_i^X,B_i^T\)，使
\[
X_c\overset d=\sum_i B_i^X E_i,\qquad
T_c\overset d=\sum_i B_i^T E_i,
\]
其中
\[
P(B_i^X=1)=p_i^X=1-\frac{\alpha_i}{\rho_i},\qquad
P(B_i^T=1)=p_i^T=1-\frac{\alpha_i}{\sigma_i}.
\]
由于 \(p_i^T\ge p_i^X\)，可用共同均匀随机数耦合成
\[
B_i^X\le B_i^T\quad\text{逐项成立},
\]
从而 \(X_c\le T_c\) 几乎处处。这是一个**同速率、不同激活概率**的 Bernoulli--指数和表示。它与 \(T_c\overset d=X_c+Z_c\) 的独立卷积表示等价，但两种表示的耦合结构不同：这里的差值与 \(X_c\) 不独立。

进一步令
\[
\tau_i^X=-\frac{\log p_i^X}{\alpha_i},\qquad
\tau_i^T=-\frac{\log p_i^T}{\alpha_i}.
\]
则
\[
B_i^X E_i\overset d=(E_i-\tau_i^X)_+,\qquad
B_i^T E_i\overset d=(E_i-\tau_i^T)_+,
\]
并且
\[
\tau_i^T\le\tau_i^X.
\]
所以两个 CDF 可统一写成
\[
F(t;\boldsymbol\tau)
=P\!\left(\sum_i(E_i-\tau_i)_+\le t\right),
\qquad
F_X(t)=F(t;\boldsymbol\tau^X),\quad
F_T(t)=F(t;\boldsymbol\tau^T).
\]
这将 \(F_T/F_X\) 的比较转化为同一个凸次水平集概率对阈值向量的比较。

### 5.5 相对 CDF 曲率的参数路径归约

对任意固定阶段 \(i\)，令 \(F_{-i}(t)\) 为去掉该阶段后的和的 CDF，并令
\[
H_i(t)=\int_0^t \alpha_i e^{-\alpha_i(t-u)}F_{-i}(u)\,du.
\]
则精确地
\[
F(t;\boldsymbol\tau)
=(1-p_i)F_{-i}(t)+p_iH_i(t),
\qquad p_i=e^{-\alpha_i\tau_i}.
\]
因此
\[
\boxed{
\partial_{\tau_i}\log F
=\alpha_i p_i
\frac{F_{-i}-H_i}
{(1-p_i)F_{-i}+p_iH_i}
\ge0.
}
\]
沿着阈值线段
\[
\boldsymbol\tau(s)
=\boldsymbol\tau^X+s(\boldsymbol\tau^T-\boldsymbol\tau^X),
\qquad 0\le s\le1,
\]
有恒等式
\[
\boxed{
\left(\log\frac{F_T}{F_X}\right)''
=
\int_0^1\sum_i
(\tau_i^T-\tau_i^X)
\,\partial_{tt\tau_i}\log F(t;\boldsymbol\tau(s))\,ds.
}
\]
因 \(\tau_i^T-\tau_i^X\le0\)，若能证明相应谱阈值路径上
\[
\partial_{tt\tau_i}\log F(t;\boldsymbol\tau)\ge0,
\]
就得到
\[
h_T'(t)-h_X'(t)
=\left(\log\frac{F_T}{F_X}\right)''\le0.
\]
这会给出一个比单独使用 \(F_X,F_T\) 对数凹性更强的比较：不仅 \(h_T\ge h_X\)，而且 \(h_T-h_X\) 非增。

**目前没有把上述混合三阶导数符号写成定理。** Prékopa 定理只给出 \(\log F\) 对 \((t,\boldsymbol\tau)\) 的联合凹性，不能单独推出 \(\partial_{tt\tau_i}\log F\ge0\)。所需符号必须利用谱阈值 \(\tau_i=-\alpha_i^{-1}\log(1-\alpha_i/\rho_i)\) 的额外结构证明。即使得到 \(h_T'\le h_X'\)，还需要控制
\[
\partial_t^2\log(h_T-\alpha)
-\partial_t^2\log(h_X+\beta)
\]
才能关闭完整的 \((\log r)''<0\)。

## 6. 真正剩下的 Robin 专属曲率不等式

现在正值区间已确定为 \((0,t_*)\)。在该区间，
\[
r(t)=C\frac{F_T(t)}{F_X(t)}
\frac{h_T(t)-\alpha}{h_X(t)+\beta},
\qquad C=\frac baM_Z>0.
\]
因此对数曲率可写成更直接的形式：
\[
\boxed{
\begin{aligned}
(\log r)''
={}&h_T'-h_X'
+\frac{h_T''}{h_T-\alpha}
-\frac{(h_T')^2}{(h_T-\alpha)^2}\\
&-\frac{h_X''}{h_X+\beta}
+\frac{(h_X')^2}{(h_X+\beta)^2}.
\end{aligned}
}
\]
需要证明该式在 \(0<t<t_*\) 严格小于零，或找到足以推出 \(r\) 单峰的较弱不等式。

目前已知：
\[
h_T'\le0,\qquad h_X'\le0,\qquad h_T\ge h_X.
\]
但这些一阶关系还不能决定 \(h_T'-h_X'\) 或两个二阶曲率组合的总符号。特别是，不能从 CDF 对数凹直接推出 \(h_X''\ge0\)，也不能把任意正卷积误认为保留对数凹性。

所以，当前缺口已经从两个条件缩小为一个明确命题：
\[
\boxed{
\mathcal C(t):=(\log r(t))''<0,
\qquad 0<t<t_*.
}
\]
若该不等式成立，则 \(r\) 在唯一正值区间严格对数凹；结合第 4 节的传递引理，立即得到 \(Q\) 至多单峰，进而 \(J\) 至多有两个正时间零点。

## 7. 文献位置与适用边界

- Bondesson, *Classes of Infinitely Divisible Distributions and Densities* (1981), [DOI:10.1007/BF00533713](https://doi.org/10.1007/BF00533713)：给出指数混合卷积类及其与首次通过时间分布的联系；GGC 子类具有单峰性，但这并不直接推出当前两个响应的比值对数凹。
- Keilson, *On the Unimodality of Passage Time Densities in Birth-Death Processes* (1981), [DOI:10.1111/j.1467-9574.1981.tb00710.x](https://doi.org/10.1111/j.1467-9574.1981.tb00710.x)：首次通过时间密度的单峰性及其卷积分解提供可借鉴结构，但没有给出本模型所需的导数比值定理。
- Di Crescenzo and Ricciardi, *Comparing First-Passage Times for Semi-Markov Skip-Free Processes* (1996), [DOI:10.1016/0167-7152(95)00227-8](https://doi.org/10.1016/0167-7152(95)00227-8)：比较半马尔可夫 skip-free 过程的首次通过时间；并未给出分层 Robin 扩散的导数比值对数凹性。
- Sangüesa, Badía and Cha, *Preservation of ageing classes in deterioration models with independent increments* (2014), [DOI:10.1186/1029-242x-2014-200](https://doi.org/10.1186/1029-242x-2014-200)：研究独立正增量过程下 CDF 对数凹（DRHR）等性质的保持；本文对倾斜阶段和的 CDF 直接使用 Prékopa 定理，不依赖其过程假设。
- Xu and Balakrishnan, *On the Convolution of Heterogeneous Bernoulli Random Variables* (2011), [DOI:10.1017/S0021900200008391](https://doi.org/10.1017/S0021900200008391)：研究异质 Bernoulli 和的似然比序与反向风险率序。这里的倾斜谱变量是带指数权重的 Bernoulli 和，故该文提供离散激活变量比较的相关方法；其结论并不直接覆盖连续加权和或本节所需的混合三阶导数符号。
- Yu, *Stochastic Ordering of Exponential Family Distributions and Their Mixtures* (2009), [DOI:10.1239/jap/1238592127](https://doi.org/10.1239/jap/1238592127)：以相对对数凹序比较指数族及其混合，并处理 Gamma 卷积；可作为后续比较 (F_T/F_X) 的工具背景，但没有直接给出当前 Robin 谱乘积所需的曲率不等式。

## 8. 当前结论

1. 两个导数响应的精确半谱卷积恒等式已经得到。
2. \(G_L'>0\)、\(G_L>0\)、\(G_R>0\) 已由正等待表示直接证明。
3. \(r\) 的正值集合是唯一的初始区间，且 \(h_T\ge h_X\) 已由对数凹核推出。
4. “\(r\) 单峰 \(\Rightarrow Q\) 单峰 \(\Rightarrow J\) 至多两个零点”的传递环节已经严格证明。
5. 尚未证明的是 Robin 特定的 \(\mathcal C(t)<0\)，或其他足以保证 \(r\) 单峰的曲率比较。

本文件不把尚未证明的比值单峰性写成定理，也不改变已有的静态打平局部结论。下一步只研究第 6 节的 Robin 专属曲率不等式，不再重复已有的全时间局部证书。
