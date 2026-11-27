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

## 5. 尚缺的 Robin 专属不等式

令
\[
H(t)=(\nu*k)(t),\qquad
u(t)=(\log H(t))',\qquad
v(t)=(\log k(t))'.
\]
在 \(r(t)>0\) 的区间内，
\[
r(t)=\frac ba\,\frac{H(t)}{k(t)}
\frac{u(t)+a}{v(t)+b},
\]
并且 \(v(t)+b>0\)，因为 \(k\) 的严格对数凹性给出 \(v'(t)<0\)，而 \(v(t)\to-c\) 且 \(b>c\)。

对数曲率的精确表达式为
\[
\boxed{
(\log r)''
=u'-v'
+\frac{u''}{u+a}
-\frac{(u')^2}{(u+a)^2}
-\frac{v''}{v+b}
+\frac{(v')^2}{(v+b)^2}.
}
\]
现有结果只给出 \(k\) 的严格对数凹性，即 \(v'<0\)。它并不自动给出 \(H=\nu*k\) 的对数凹性，也不能仅凭 \(\nu\) 是正测度就断言 \(H/k\) 对数凹。任意正测度卷积并不普遍保持对数凹性。

因此，不能把“\(k\) 对数凹”直接替换成“\(r\) 单峰”。需要利用 \(\nu\) 的特殊 Robin 谱结构，证明上式严格小于零，并证明 \(u+a\) 只变号一次。一般正常半谱排序不足以自动提供这两个结论；已有静态打平风险率反例也说明，不能仅凭 \(L_i\ge R_i\) 断言任意相关风险率只过阈值一次。

## 6. 文献位置与适用边界

- Bondesson, *Classes of Infinitely Divisible Distributions and Densities* (1981), [DOI:10.1007/BF00533713](https://doi.org/10.1007/BF00533713)：给出指数混合卷积类及其与首次通过时间分布的联系；GGC 子类具有单峰性，但这并不直接推出当前两个响应的比值对数凹。
- Keilson, *On the Unimodality of Passage Time Densities in Birth-Death Processes* (1981), [DOI:10.1111/j.1467-9574.1981.tb00710.x](https://doi.org/10.1111/j.1467-9574.1981.tb00710.x)：首次通过时间密度的单峰性及其卷积分解提供可借鉴结构，但没有给出本模型所需的导数比值定理。
- Di Crescenzo and Ricciardi, *Comparing First-Passage Times for Semi-Markov Skip-Free Processes* (1996), [DOI:10.1016/0167-7152(95)00227-8](https://doi.org/10.1016/0167-7152(95)00227-8)：比较半马尔可夫 skip-free 过程的首次通过时间；并未给出分层 Robin 扩散的导数比值对数凹性。

## 7. 当前结论

1. 两个导数响应的精确半谱卷积恒等式已经得到。
2. \(G_L'>0\)、\(G_L>0\)、\(G_R>0\) 已由正等待表示直接证明。
3. “\(r\) 单峰 \(\Rightarrow Q\) 单峰 \(\Rightarrow J\) 至多两个零点”的传递环节已经严格证明。
4. 尚未证明的是 Robin 特定的 \(r\) 正值区间单连通性和严格对数凹性。

本文件不把尚未证明的比值单峰性写成定理，也不改变已有的静态打平局部结论。下一步只研究第 5 节的对数曲率不等式，不再重复已有的全时间局部证书。
