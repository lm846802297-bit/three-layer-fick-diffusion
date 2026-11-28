# Robin 曲率短时间区间：反射路径的显式余项证书

## 1. 目标与基准

本证书处理静态等阻力三层 Robin 基准
\[
w_1=\frac{25}{66},\quad D_1=25,\qquad
w_2=\frac4{33},\quad D_2=\frac14,
\]
右半区长度 \(1/2\)、扩散率 \(1\)，Robin 参数 \(\kappa=4\)。令
\[
\delta=\ell_R-\ell_L=\frac2{11},\qquad
\ell_L=\frac7{22},\quad \ell_R=\frac12,
\]
以及 \(c=\lambda_2\)、\(q=\sqrt{s-c}\)。已有区间根证书给出
\[
L_1<14.24,\qquad c<15.81,\qquad 0<R_1-L_1<6.73.
\tag{1.1}
\]
本节证明
\[
\boxed{(\log f_Z)''(t)<0,\qquad 0<t\leq t_0:=0.0015.}
\tag{1.2}
\]

## 2. 精确传递矩阵分解

令
\[
a=\frac5{66},\qquad b=\frac{16}{66},\qquad a+b=\ell_L.
\]
对 \(k=iq\)，左半区的精确特征函数为
\[
\begin{aligned}
Y_L(iq)={}&\cosh(aq)\cosh(bq)+10\sinh(aq)\sinh(bq)\\
&+\frac{\kappa}{q}\left[\frac15\cosh(bq)\sinh(aq)
+2\sinh(bq)\cosh(aq)\right].
\end{aligned}
\]
直接整理得
\[
Y_L(iq)=\frac{11}{4}e^{(a+b)q}
\left(1+\frac{\kappa}{5q}\right)N_L(q),
\]
其中
\[
N_L(q)=1-\frac9{11}\frac{q-\kappa/5}{q+\kappa/5}e^{-2aq}
-\frac9{11}e^{-2bq}
+\frac{q-\kappa/5}{q+\kappa/5}e^{-2(a+b)q}.
\tag{2.1}
\]
右半区有
\[
Y_R(iq)=\frac12e^{q/2}\left(1+\frac\kappa q\right)N_R(q),
\]
\[
N_R(q)=1+\frac{q-\kappa}{q+\kappa}e^{-q}.
\tag{2.2}
\]

由于首个左右半谱根不同，去掉首根的规范乘积贡献因子 \(L_1/R_1\)。因此延迟变量的精确 Laplace 变换为
\[
\boxed{
\Phi_Z(s)=K e^{-\delta q}H(q)\frac{N_L(q)}{N_R(q)},
}
\tag{2.3}
\]
其中
\[
K=\frac{11C_0L_1}{2R_1}>0,
\qquad
H(q)=\frac{q+\kappa/5}{q+\kappa}
\frac{q^2+R_1}{q^2+L_1}.
\tag{2.4}
\]
这里保留了 \(L_1/R_1\)；不能把两个半谱首根当作共同根。

将主项与反射项分开：
\[
\Phi_0(s)=K e^{-\delta q}H(q),\qquad
\Phi_{\rm ref}(s)=K e^{-\delta q}H(q)
\left(\frac{N_L(q)}{N_R(q)}-1\right).
\tag{2.5}
\]

## 3. 主项的三个正积分与导数控制

记
\[
g_\delta(t)=\frac{\delta}{2\sqrt\pi}t^{-3/2}
 e^{-\delta^2/(4t)},
\qquad
h_\kappa=\mathcal L^{-1}\left\{\frac1{\sqrt p+\kappa}\right\}.
\]
由
\[
\frac1{\sqrt p+\kappa}=\int_0^\infty e^{-\kappa u}e^{-u\sqrt p}\,du,
\qquad
\mathcal L^{-1}\{e^{-u\sqrt p}\}=g_u(t),
\]
定义
\[
r(t)=\frac{g_\delta*h_\kappa}{g_\delta},\qquad
\ell(t)=\frac{g_\delta*e^{-L_1\cdot}}{g_\delta},\qquad
j(t)=\frac{g_\delta*h_\kappa*e^{-L_1\cdot}}{g_\delta}.
\]
于是
\[
\boxed{f_0(t)=K e^{ct}g_\delta(t)H_0(t),\qquad
H_0=1-a_0r+b_0\ell-a_0b_0j,}
\tag{3.1}
\]
其中 \(a_0=4\kappa/5=3.2\)、\(b_0=R_1-L_1<6.73\)。

### 3.1 单卷积 \(r\)

作变量替换 \(x=((\delta+u)^2-\delta^2)/(4t)\)，得到精确积分
\[
r(t)=\frac{2t}{\delta}\int_0^\infty
 e^{-x-\kappa(\sqrt{\delta^2+4tx}-\delta)}\,dx.
\tag{3.2}
\]
令 \(d=\sqrt{\delta^2+4tx}\)。有
\[
r\leq\frac{2t}{\delta},\qquad
|r'|\leq\frac2\delta+\frac{4\kappa t}{\delta^2},
\]
以及
\[
|r''|\leq\frac{8\kappa}{\delta^2}
+\frac{16\kappa^2t}{\delta^3}
+\frac{16\kappa t}{\delta^4}.
\]
因为
\[
\frac{2\kappa t_0}{\delta}+\frac{2t_0}{\delta^2}<1,
\]
故在 \(0<t\leq t_0\) 上
\[
\boxed{r\leq\frac{2t}{\delta},\quad |r'|<12,
\quad |r''|\leq\frac{16\kappa}{\delta^2}.}
\tag{3.3}
\]

### 3.2 单卷积 \(\ell\)

作变量替换 \(x=u/(t-u)\)，再令 \(y=\delta^2x/(4t)\)，得到
\[
\ell(t)=\frac{4t^2}{\delta^2}\int_0^\infty e^{-y}F_\ell(t,y)\,dy,
\]
\[
F_\ell(t,y)=\left(1+\frac{4ty}{\delta^2}\right)^{-1/2}
\exp\left[-\frac{4L_1t^2y}{\delta^2+4ty}\right].
\]
有 \(0<F_\ell\leq1\)，并且
\[
|\partial_t\log F_\ell|\leq\frac1{2t}+L_1,
\qquad
|\partial_t^2\log F_\ell|\leq\frac1{2t^2}+\frac{2L_1}{t}.
\]
从而
\[
|\ell'|\leq\frac{10t}{\delta^2}
+\frac{4L_1t^2}{\delta^2},
\]
\[
|\ell''|\leq\frac{19}{\delta^2}
+\frac{28L_1t}{\delta^2}
+\frac{4L_1^2t^2}{\delta^2}.
\]
利用 (1.1)，在整个短时区间上
\[
\boxed{|\ell'|<1,\qquad
|\ell''|\leq\frac{\delta^2/2+4t_0}{t^2}.}
\tag{3.4}
\]

### 3.3 双卷积 \(j\)

从 \(g_\delta*h_\kappa=\int_0^\infty e^{-\kappa u}g_{\delta+u}\,du\) 出发，依次令
\[
y=\frac{(\delta+u)^2-\delta^2}{4t},
\qquad
z=\frac{(\delta+u)^2x}{4t},
\]
得到
\[
j(t)=\frac{8t^3}{\delta}\int_0^\infty\!\int_0^\infty
 e^{-y-z}F_j(t,y,z)\,dy\,dz,
\tag{3.5}
\]
其中
\[
\begin{aligned}
F_j={}&e^{-\kappa(d-\delta)}d^{-2}
\left(1+\frac{4tz}{d^2}\right)^{-1/2}\\
&\times\exp\left[-\frac{4L_1t^2z}{d^2+4tz}\right],
\qquad d^2=\delta^2+4ty.
\end{aligned}
\]
故 \(0<F_j\leq\delta^{-2}\)，并直接得到
\[
0\leq j\leq\frac{8t^3}{\delta^3}
\leq\frac{16t^3}{\delta^3}.
\tag{3.6}
\]

令
\[
A_y=\frac{2\kappa}{\delta}+\frac4{\delta^2},
\qquad
A_z=\frac2{\delta^2}+\frac{8L_1t_0}{\delta^2}.
\]
对 \(\phi_j=\log F_j\) 求导，使用 \(d\geq\delta\)，有
\[
|\partial_t\phi_j|\leq A_yy+A_zz.
\]
再用两个独立单位指数变量的矩
\(\mathbb E Y=\mathbb E Z=1\)、
\(\mathbb E Y^2=\mathbb E Z^2=2\)、
\(\mathbb E(YZ)=1\)，得到
\[
|j'|\leq\frac{24t^2}{\delta^3}
+\frac{8t^3(A_y+A_z)}{\delta^3}
\leq\frac{2t}{\delta},
\tag{3.7}
\]
其中最后一步的无量纲比值在 \(t=t_0\) 处小于 \(0.608\)。

对二阶导数，\(\phi_j\) 的二阶导数可由
\[
|\partial_t^2\phi_j|
\leq \frac{4\kappa y^2}{\delta^3}
+\frac{16y^2+16yz+8z^2}{\delta^4}
+\frac{8L_1z}{\delta^2}
\]
控制。将其与 \((\partial_t\phi_j)^2\) 一起积分，得到
\[
|j''(t)|\leq
\frac{48t}{\delta^3}
+\frac{48t^2(A_y+A_z)}{\delta^3}
+\frac{8t^3}{\delta^3}
\left[2A_y^2+2A_z^2+2A_yA_z
+\frac{8\kappa}{\delta^3}+\frac{80}{\delta^4}
+\frac{8L_1}{\delta^2}\right].
\]
右侧在 \(t\leq t_0\) 上递增，其在 \(t_0\) 的上界为
\[
\boxed{|j''(t)|<16.873
<12+\frac{\delta}{t_0}+\frac8\delta.}
\tag{3.8}
\]

### 3.4 主项的统一曲率界

将 (3.3)、(3.4)、(3.6)–(3.8) 代入 (3.1)，得到
\[
H_0(t)>0.9470065,
\qquad |H_0'(t)|<45.486,
\]
以及
\[
t_0^2|H_0''(t)|<0.174146,
\qquad
 t_0^2| (\log H_0)''(t)|<0.189082<0.20.
\]
因此
\[
| (\log H_0)''(t)|<\frac{0.20}{t^2}.
\]
因为
\[
(\log g_\delta)''=-\frac{\delta^2}{2t^3}
+\frac{3}{2t^2},
\]
并且
\[
\frac{\delta^2}{2}-1.7t_0>0.0139,
\]
所以主项满足
\[
\boxed{(\log f_0)''(t)<-\frac{0.0139}{t^3},
\qquad 0<t\leq t_0.}
\tag{3.9}
\]
另外还得到后续商估计所需的导数界。首先，
\[
t^2(\log f_0)'=\frac{\delta^2}{4}-\frac32t
+t^2\left(c+\frac{H_0'}{H_0}\right).
\]
利用 \(c<15.81\)、\(|H_0'/H_0|<45.486/0.9470065\)，其下界在 \(t=t_0\) 仍严格为正；其上界小于 \(0.009\)。故
\[
0<(\log f_0)'(t)<\frac{0.009}{t^2}.
\]
其次，(3.9) 已经证明 \((\log f_0)''<0\)，而
\[
\frac{3}{2t^2}+(\log H_0)''>\frac{1.3}{t^2}>0.
\]
因此曲率的绝对值小于其负的主奇异项绝对值：
\[
\boxed{|(\log f_0)'(t)|<\frac{0.009}{t^2},
\qquad |(\log f_0)''(t)|<\frac{\delta^2}{2t^3}<\frac{0.019}{t^3}.}

\]
\tag{3.10}

## 4. 反射项的 Bromwich 统一界

令
\[
d_0=\delta+2a=\frac13,
\qquad
\alpha=\frac{d_0^2-\delta^2}{4}
=0.0195133149678604\ldots.
\tag{4.1}
\]
在 \(\Re q>0\) 内，
\[
\left|\frac{q-\kappa/5}{q+\kappa/5}\right|\leq1,
\qquad
\left|\frac{q-\kappa}{q+\kappa}\right|\leq1.
\]
由 (2.1)–(2.2)，
\[
|N_L-1|\leq\frac{29}{11}e^{-2a\Re q},
\qquad |N_R-1|\leq e^{-\Re q}.
\]
取 \(t\leq t_0\) 的 Bromwich 直线
\[
\Re p=\sigma=\frac{d_0^2}{4t^2},\qquad q=\sqrt p.
\]
则
\[
\Re q\geq\sqrt\sigma\geq\frac{d_0}{2t_0}>111.
\]
因此
\[
|e^{2aq}(N_L-N_R)|<3.
\]
同时，\(H\) 的第一个有理因子模长不超过 1，第二个满足
\[
\left|\frac{q^2+R_1}{q^2+L_1}\right|
\leq1+\frac{R_1-L_1}{\sigma+L_1}
<1.000546.
\]
又 \(|N_R|\geq1-e^{-\Re q}\)，故统一得到
\[
\boxed{
\left|H(q)e^{2aq}\left(\frac{N_L}{N_R}-1\right)\right|
<3.001.
}
\tag{4.2}
\]

### 4.1 一个显式的 Bromwich 积分引理

定义
\[
J_n(A)=\int_{-\infty}^{\infty}
(1+v^2)^{n/2}
 e^{-A(\Re\sqrt{1+iv}-1)}\,dv,
\qquad n=0,1,2.
\]
令 \(w(v)=\Re\sqrt{1+iv}\)。由
\[
w(v)-1=
\frac{v^2}{2(\sqrt{1+v^2}+1)(w(v)+1)},
\]
以及 \(\sqrt2<1.415\)，可得
\[
w(v)-1\geq\begin{cases}
v^2/11,&|v|\leq1,\\
\sqrt{|v|}/11,&|v|\geq1.
\end{cases}
\]
因此
\[
J_n(A)\leq
\sqrt{\frac{11\pi}{A}}
+\frac{\sqrt\pi}{2}\left(\frac{11}{A}\right)^{3/2}
+960\left(\frac{11}{A}\right)^6.
\]
对 \(A\geq37\)，右侧小于 \(64/\sqrt A\)。所以
\[
\boxed{J_n(A)<\frac{64}{\sqrt A},
\quad n=0,1,2,\quad A\geq37.}
\tag{4.3}
\]

### 4.2 反射密度及其导数

在 (4.2) 中取 \(d_0=1/3\)，并在 Bromwich 线上令
\[
p=\sigma(1+iv),\qquad
\sigma=\frac{d_0^2}{4t^2},\qquad
A=\frac{d_0^2}{2t}.
\]
则
\[
\sigma t-d_0\Re\sqrt p
=-\frac{d_0^2}{4t}-A(\Re\sqrt{1+iv}-1).
\]
对 \(n\) 阶时间导数，Laplace 积分中出现 \((p+c)^n\)，且
\[
|p+c|^n\leq
\left(1+\frac{c}{\sigma}\right)^n|p|^n.
\]
在 \(t\leq t_0\) 上，
\[
1+\frac{c}{\sigma}<1.002,
\qquad e^{ct}<e^{15.81t_0}.
\]
结合 (4.2)–(4.3)，并用 \(t_0<\alpha/4\) 保证各指数幂因子在 \((0,t_0]\) 上递增，得到
\[
\frac{|f_{\rm ref}^{(n)}(t)|}{K g_\delta(t)}
\leq C_n t^{-2n}e^{-\alpha/t},
\qquad n=0,1,2,
\tag{4.4}
\]
其中可取
\[
C_0<73,\qquad C_1<2.1,\qquad C_2<0.057.
\tag{4.5}
\]

## 5. 对数曲率的最终误差比较

写
\[
f_Z=f_0(1+\varepsilon),
\qquad \varepsilon=\frac{f_{\rm ref}}{f_0}.
\]
由 \(H_0>0.94\)、(3.10)、(4.4)–(4.5)，以及
\(t^{-m}e^{-\alpha/t}\) 对 \(m\leq7\) 在 \((0,t_0]\) 上递增，得到
\[
\sup|\varepsilon|<1.73985\times10^{-4},
\]
\[
\sup|\varepsilon'|<2.920405,
\qquad
\sup|\varepsilon''|<4.83938\times10^4.
\]
因此
\[
\begin{aligned}
\left|(\log(1+\varepsilon))''\right|
&\leq\frac{\sup|\varepsilon''|}{1-\sup|\varepsilon|}
+\frac{(\sup|\varepsilon'|)^2}{(1-\sup|\varepsilon|)^2}\\
&<4.8411\times10^4.
\end{aligned}
\tag{5.1}
\]
而 (3.9) 在整个短时区间上给出统一上界
\[
(\log f_0)''(t)<-4.1185185\times10^6.
\]
故
\[
\boxed{
(\log f_Z)''(t)<-4.0701\times10^6<0,
\qquad 0<t\leq0.0015.
}
\tag{5.2}
\]

## 6. 与已认证时间区间的连接

此前已归档的区间证书分别给出：

- \([0.0015,0.004]\)：\(N=30\)，根证书延伸至 \(R_{70}\)，并对 \(30<j\leq70\) 逐项控制、对 \(j>70\) 使用 Gaussian 谱尾；
- \([0.004,0.010]\)：有限谱与 Gaussian 尾的严格负性；
- \([0.010,0.012]\)：\(Q_Z'>0\)，唯一零点括在 \((0.01102,0.01103)\)；
- \([0.012,\infty)\)：首对谱项支配，严格正性。

因此短时证书 (5.2) 正好关闭此前唯一缺失的 \((0,0.0015)\) 区间。对当前 \(\kappa=4\) 静态打平基准，完整时间轴的曲率符号已全部覆盖。

## 7. 可复核文件与运行环境

- robin_short_time_reflection_certificate.py：70 位区间算术常数验收。
- robin_short_time_reflection_certificate_input.json：精确基准、根括区与验收参数。
- robin_short_time_reflection_certificate_output.txt：本机原生运行输出。
- robin_short_time_reflection_certificate_provenance.json：源码、输入、输出与本证明文档的 SHA-256。

运行环境：Python 3.13.5、mpmath 1.3.0、Linux x86_64。短时证明不重新计算全谱；仅使用此前已认证的 \(L_1,R_1,c\) 根区间及其导出的粗上界。
