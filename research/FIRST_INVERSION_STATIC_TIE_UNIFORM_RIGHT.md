# 静态打平：用均匀右半区的精确谱封住风险率的晚时段

继续研究[时间域单零点判据](FIRST_INVERSION_STATIC_TIE_TIME.md)。记 \(a=L_1<b=R_1\)，对 \(i\ge2\) 有 \(L_i\ge R_i\)。令 \(Z\) 是这些正常谱对的等待和，\(S=\Pr(Z>t)\)，\(f=-S'\)，\(h=f/S\)。本笔记证明 \(h>a\) 的显式晚时界，**不把此界当作 \(J\) 的零点上限**。

## 1. 剥离最慢的有效正常模

取第一个严格正常的指标 \(j=\min\{i\ge2:L_i>R_i\}\)，因光学长度严格有差，该指标存在。记 \(m=R_j,\ q=1-m/L_j>0\)，将 \(Z=Z_j+T\) 分解；\(Z_j\) 的非零概率为 \(q\)，非零时呈速率 \(m\) 的指数律，\(T\) 为后续独立等待和。记 \(M_\eta=\mathbb E e^{\eta T}\)。给定 \(T\) 特别保留仍未进入第 \(j\) 段的概率，得到
\[
 S(t)=S_T(t)+q e^{-mt}\mathbb E[e^{mT};T\le t],\quad
 f(t)=(1-q)f_T(t)+qm e^{-mt}\mathbb E[e^{mT};T\le t].
\]
故对 \(m<\gamma<R_{j+1}\) 有
\[
 f-aS\ \ge\ q(m-a)M_m e^{-mt}
  -[a+q(m-a)]M_\gamma e^{-\gamma t}.
\]
这里分别用 \(S_T\le M_\gamma e^{-\gamma t}\) 和
\(\mathbb E[e^{mT};T>t]\le M_\gamma e^{-(\gamma-m)t}\)，且只丢弃了非负的 \((1-q)f_T\)。严格光学长度差也使移除有限多个因子后的等待和无原子；密度等式可先按测度推导，再用下述传递矩阵反演的光滑性写成逐点等式。

## 2. 均匀 Robin 右半区将无穷尾积变成初等函数

右半长度 \(1/2\)、扩散率 1，右特征值 \(R_i=k_i^2\) 的精确方程是
\[
 k_i=(2i-1)\pi+2\arctan(\kappa/k_i),\qquad
 (2i-1)\pi<k_i<2i\pi .
\]
用解 \(y(0)=1,\ y'(0)=\kappa\) 的中点迹，归一化整函数的精确乘积为
\[
 P_R(-\gamma)=\frac{\cos(\sqrt\gamma/2)+
      (\kappa/\sqrt\gamma)\sin(\sqrt\gamma/2)}{1+\kappa/2}
   =\prod_{i\ge1}(1-\gamma/R_i).
\]
后续所有 \(1-\gamma/L_i\) 不超过 1，故
\[
 M_\gamma \le
 U_\gamma:=\frac{\prod_{i=1}^j(1-\gamma/R_i)}{P_R(-\gamma)}.
\]
分子与分母有同号，商严格正。无须近似高频根，也无须知道任何 \(i>j\) 的左半谱。若已另有几个后续 \(L_i\) 的认证区间，可乘上对应的 \(1-\gamma/L_i\) 来收紧；若选用 \(\gamma=(2j+1/2)^2\pi^2\)，上面的右根间隔保证其永远落在目标谱隙内。由 \(M_m\ge1\) 得到实用的严格界：
\[
 \boxed{\qquad h(t)>a
 \quad\text{只要}\quad
 t>\frac{\log\{U_\gamma[1+a/(q(m-a))]\}}{\gamma-m} .\qquad}
\]
右边依赖慢模的实际概率 \(q\)；当 \(q\to0\) 时本估计并不一致，不能当成对整张参数面的统一截止。

## 3. 早时端的解析符号

令 \(\delta=\ell_R-\ell_L>0\) 为两半光学长度差；静态打平时严格 Cauchy–Schwarz 给出该正性。左半两层从壁起扩散率 \(D_{\rm w},D_{\rm m}\)，固定正厚度、正扩散率。对 \(\operatorname{Re}s>0\) 的大 \(s\) 用传递矩阵展开，取 \(\operatorname{Re}\sqrt s>0\)，有
\[
 \widehat\nu(s)=C e^{-\delta\sqrt{s}}
     \bigl[1+O(|s|^{-1/2})+O(e^{-\eta\operatorname{Re}\sqrt s})\bigr],
 \quad C=\frac{a}{2b}\bigl(1+\sqrt{D_{\rm w}/D_{\rm m}}\bigr)>0 .
\]
系数来自左半前向项 \((1+\sqrt{D_{\rm w}/D_{\rm m}})/4\)、右半前向项 \(1/2\) 和移去首谱因子的高频常数 \(a/b\)。这也给密度及任意固定阶导数的 Bromwich 绝对收敛。与[前次短时反演](FIRST_INVERSION_PHYSICAL_INPUT_NO_STACKING.md)相同，取 \(\operatorname{Re}s=\delta^2/(4t^2)\)；在该线上对余项绝对积分（\(\operatorname{Re}\sqrt{1+iv}-1\gtrsim\min(v^2,\sqrt{|v|})\)）得到
\[
 f(t)=\frac{C\delta}{2\sqrt\pi\,t^{3/2}}e^{-\delta^2/(4t)}
 (1+O(t)), \qquad t\downarrow0.
\]
乘 \(s\) 后作同样的估计，领先项是上述首项对 \(t\) 的导数，正且相对误差仍为 \(O(t)\)。因此 \(h\to0\) 且 \(h'>0\) 在充分小的正时间成立。有限段因子乘除只改变高频常数和代数误差，本论证同样适用于 \(T\) 去掉的有限前缀。这里没有为任意参数算出早时端的显式下截止。

## 4. 严格算例及剩余条件

从左壁到中点依次是厚 \(25/66\) 的 \(D=25\) 快层和厚 \(4/33\) 的 \(D=1/4\) 慢层；右半 \(D=1\)，两壁 \(\kappa=4\)。左、右单位横截面扩散阻力同为 \(1/2\)，这是有理恒等式。以区间三角函数认证前 3 个右根处的左连续相位；由光学间隙，\(i\ge4\) 自动正常。因此确实只有首模逆序。取 \(\gamma=183\)，区间证书得到 \(a\simeq14.2331,\ m\simeq103.5097,\ q\simeq0.41106,\ U_\gamma\simeq17.40665\)，右边界约 \(0.04006281\)。于是对 **全部 \(t\ge0.041\)**，\(h(t)>a\)，并非对许多时间点逐个采样。见配套 Python 和原始输出。

结合早时 \(h<a\)，所有再次穿越阈值的可能性都集中在两端之间。尚未证明该中间窗只有一次阈值上穿，尚未证明整个打平面排除三个 \(J\) 零点，也没有得到真实 \(\Delta\) 反超次数的全域界。
