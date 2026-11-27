# 静态打平面：用两侧精确特征函数封闭正常谱尾矩

本笔记是 [均匀右半区风险率晚时界](FIRST_INVERSION_STATIC_TIE_UNIFORM_RIGHT.md) 的严格加强。此前只用右半区特征函数给
\(M_\gamma\) 上界，并以 \(M_m\ge1\) 处理慢模矩。本轮同时使用左侧两层和右侧均匀半区的精确特征函数，得到 \(M_\gamma\) 与 \(M_m\) 的无截断表达式。

**基准结论：**
\[
\boxed{\ h(t)>a=L_1\quad\text{对所有 }t\ge 0.0171.\ }
\]
原来的严格晚时截止为 \(t\ge0.041\)。新截止来自完整半谱乘积，不使用高阶本征值截断。它仍然没有关闭 \(0<t<0.0171\) 的中间时间窗，也不推出 \(J\) 的全域唯一零点。

## 1. 正常谱等待和原有下界

假设仅首模逆序，\(a=L_1<b=R_1\)，且 \(L_i\ge R_i\) 对 \(i\ge2\)。令
\[
m=R_2,\qquad q=1-\frac{m}{L_2}>0,\qquad Z=Z_2+T,
\]
其中 \(Z_2=0\) 的概率为 \(1-q\)，以概率 \(q\) 服从速率 \(m\) 的指数分布；\(T\) 是 \(i\ge3\) 的独立正常谱等待和。记
\[
S(t)=\Pr(Z>t),\qquad f(t)=-S'(t),\qquad
M_\eta=\mathbb E e^{\eta T}.
\]
对 \(m<\gamma<R_3\)，有
\[
M_\eta=\prod_{i\ge3}\frac{1-\eta/L_i}{1-\eta/R_i}
\quad(\eta=m,\gamma).
\]
分解 \(Z_2+T\)，并在正项中保留 \(T\le t\) 的部分，给出
\[
f(t)-aS(t)
\ge q(m-a)M_m e^{-mt}
-[a+q(m-a)]M_\gamma e^{-\gamma t}.
\]
因此只要
\[
t>T_\gamma:=
\frac{\log\!\left(
\displaystyle\frac{M_\gamma}{M_m}
\left[1+\frac{a}{q(m-a)}\right]\right)}
{\gamma-m},
\]
就有 \(f(t)-aS(t)>0\)，等价于 \(h(t)>a\)。此前的右侧单边估计把 \(M_\gamma\) 换成上界、把 \(M_m\) 换成 1；这里将这两个步骤都精确化。

## 2. 两侧特征函数给出的精确矩

定义归一化半区特征函数
\[
P_L(s)=\prod_{i\ge1}\left(1+\frac{s}{L_i}\right),
\qquad
P_R(s)=\prod_{i\ge1}\left(1+\frac{s}{R_i}\right).
\]
它们是相应 Robin--Dirichlet Sturm--Liouville 特征函数除以其 \(s=0\) 值后的 Hadamard 乘积。对 \(m<\gamma<R_3\)，直接剥去前两个因子：
\[
\boxed{
M_\gamma=
\frac{P_L(-\gamma)}{P_R(-\gamma)}
\frac{(1-\gamma/R_1)(1-\gamma/R_2)}
     {(1-\gamma/L_1)(1-\gamma/L_2)}.
}
\]
这个式子给出完整无穷尾积的精确值，而不只是上界。

在 \(\gamma=m=R_2\) 处，右侧特征函数有一个简单零点。利用
\[
mP_R'(-m)
=(1-m/R_1)\prod_{i\ge3}(1-m/R_i)
\]
得到
\[
\boxed{
M_m=
\frac{P_L(-m)(1-m/R_1)}
{mP_R'(-m)(1-m/L_1)(1-m/L_2)}.
}
\]
这里 \(P_R'\) 是对 \(s\) 求导。该式也正是 \(M_\gamma\) 在 \(\gamma\downarrow m\) 时的连续极限。

### 两层左半区和均匀右半区的初等表达式

左半区从壁面向中点依次经过厚度 \(w_1,w_2\)、扩散率 \(D_1,D_2\) 的两层。对 \(s=-\gamma\)，取状态向量 \((y,Dy')^\mathsf T\)，单层传递矩阵为
\[
\mathcal M_D(w;\gamma)=
\begin{pmatrix}
\cos(\sqrt{\gamma/D}\,w)&
\dfrac{\sin(\sqrt{\gamma/D}\,w)}{\sqrt{D\gamma}}\\[6pt]
-\sqrt{D\gamma}\sin(\sqrt{\gamma/D}\,w)&
\cos(\sqrt{\gamma/D}\,w)
\end{pmatrix}.
\]
令
\[
Y_L(-\gamma)=
(1,0)\,
\mathcal M_{D_2}(w_2;\gamma)
\mathcal M_{D_1}(w_1;\gamma)
\binom{1}{\kappa}.
\]
右半区均匀、长度为 \(1/2\)，所以
\[
Y_R(-\gamma)=
\cos\frac{\sqrt\gamma}{2}
+\frac{\kappa}{\sqrt\gamma}
 \sin\frac{\sqrt\gamma}{2}.
\]
本基准静态打平，且两侧零频阻力相同，因此
\[
Y_L(0)=Y_R(0)=1+\frac{\kappa}{2}=3,
\qquad
P_{L,R}(-\gamma)=Y_{L,R}(-\gamma)/3.
\]
在 \(m=R_2=k_2^2\) 处，
\[
P_R'(-m)=
-\frac{1}{2k_2}
\left.\frac{d}{dk}
\left[
\frac{\cos(k/2)+(\kappa/k)\sin(k/2)}{3}
\right]\right|_{k=k_2}.
\]
因此 \(M_\gamma\) 和 \(M_m\) 都只需有限层传递矩阵、三角函数和已认证的低模根区间。

## 3. 区间证书

输入基准为
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac{4}{33},\quad D_2=\frac14,\qquad
\kappa=4,\quad D_R=1.
\]
静态打平由精确有理恒等式给出：
\[
\frac{w_1}{D_1}+\frac{w_2}{D_2}=\frac12.
\]
左侧光学长度为
\[
\ell_L=\frac{w_1}{\sqrt{D_1}}+
       \frac{w_2}{\sqrt{D_2}}
       =\frac7{22}<\frac12=\ell_R.
\]
低模根区间沿用此前的区间隔离，并在本轮脚本中重新检查残差异号和相位排序：

| 根 | 区间 |
|---|---:|
| \(k_{L,1}\) | [3.77267640109, 3.77267640111] |
| \(k_{L,2}\) | [13.25724841026, 13.25724841029] |
| \(k_{R,1}\) | [4.57785945620, 4.57785945622] |
| \(k_{R,2}\) | [10.17397018820, 10.17397018821] |
| \(k_{R,3}\) | [16.19232720644, 16.19232720646] |

结合原检查中的高模光学相位界，只有首模逆序。取 \(\gamma=183\)，有
\[
a<m<183<R_3,\qquad q>0.
\]
50 位外向舍入区间计算给出：

| 量 | 严格区间或包络 |
|---|---:|
| \(a=L_1\) | [14.2330872273413945531881, 14.2330872274923016092321] |
| \(b=R_1\) | [20.9567972007197597184400, 20.9567972009028740966884] |
| \(m=R_2\) | [103.5096693903823434192400, 103.5096693905858228230041] |
| \(L_2\) | [175.7546354113412972732676, 175.7546354121367321778841] |
| \(R_3\) | [262.1914603604170143774736, 262.1914603610647074657316] |
| \(q=1-m/L_2\) | [0.4110558213822994746689, 0.4110558213861226805761] |
| \(M_{183}\) | [6.7409338458816155413956, 6.7409338468267028898887] |
| \(M_m\) | [2.4125847986448266153544, 2.4125847987639432302142] |
| \(T_{183}\) | [0.01704928867520349008097, 0.01704928867771593406945] |

于是 \(T_{183}<0.0171\)，严格得到
\[
\boxed{h(t)>L_1\quad\forall\,t\ge0.0171.}
\]
相对于原来的截止 \(0.041\)，新截止缩短约 \(57.4\%\)。

## 4. 文献与方法适用性

Bar-Lev 与 Letac 对自然指数族中两个成员的混合给出了风险率递增的充分条件，并讨论了 IFR 密度卷积的闭包性质：
[Increasing Hazard Rate of Mixtures for Natural Exponential Families, *Advances in Applied Probability* **44** (2012) 373--390](https://doi.org/10.1239/aap/1339878716)。
该结果主要针对两个总体的混合，不直接覆盖这里由无限多个零膨胀指数等待构成的谱等待和。更重要的是，现有 [时间域单零点分析](FIRST_INVERSION_STATIC_TIE_TIME.md) 已给出有限正常谱乘积的精确反例，说明仅凭正常谱排序不能推出风险率单调或单次阈值穿越。因此本轮没有把 IFR 闭包定理误用于本模型，而是利用两层/均匀半区的精确传递矩阵直接求出完整尾矩。

## 5. 可复现记录与结论边界

- 基础提交：5dce3e712fd4fe78d85c627c75fb82cd44b0dde4。
- 输入：本文件列出的有理几何参数、\(\gamma=183\) 和五个低模根区间；输入清单见 static_tie_uniform_right_exact_moments_input.json。
- 源码：static_tie_uniform_right_exact_moments.py；SHA-256：3a918c7c70257a3e686a47c12d1cda4ab97f574adbc36cdcfec0f570c88a8c56。
- 原生输出：static_tie_uniform_right_exact_moments_output.txt；SHA-256：65e2fa611efa84b92149d936dcc97074ee340311ce1ecf90f6f32bf123dfe3e8。
- 输入清单 SHA-256：e5fb3437c612e8fd47345f5f055318d2517101561dd14cbe82398134276513cc。
- 运行环境：Python 3.13.5，mpmath 1.3.0，Linux 6.18.57 x86_64；5 个逻辑 CPU 可见，未使用 GPU。区间精度为 50 decimal digits。
- 独立一致性检查：静态阻力和光学长度由有理数精确验证；\(M_m\) 的导数公式是 \(M_\gamma\) 在 \(\gamma\downarrow m\) 时的解析极限。高阶谱没有截断，也没有把有限模乘积近似当作证书。

**本轮完成：** 对该静态打平算例，风险率晚时正性已严格推进到 \(t\ge0.0171\)。

**仍待解决：** 早时解析结论仍为 \(h(t)<a\)；\(0<t<0.0171\) 的中间窗尚未证明只有一次上穿。因而还没有证明该算例的 \(J\) 全时间恰有一个零点，更没有证明整个静态打平面上的唯一性或实际捕获反超次数。
