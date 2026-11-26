# 双消失点的有限时间正性与文献路线

沿用 FIRST_INVERSION_WEIGHTED_LATE_PAIR.md，\(W=e^{\lambda_2t}(\partial_t+\lambda_1)\partial_t\Delta\)。
以下认证仅针对该临界点及其邻域，不声称覆盖任何其它双消失点。

## 结果

1. 在该点 \(w_2=w_3=0,w_4>0\)。145 个本征值经节点计数和区间
   传递矩阵确认，\([0.002,0.1]\) 的区间 Taylor 界证明
   \(e^{(\lambda_4-\lambda_2)t} W(t)>0\)；对 \(t\ge0.1\) 用第一非零项和
   统一尾界证明同一结论。没有直接认证 \(0<t<0.002\) 的 \(W\)。
2. 设 \(\alpha=\lambda_3-\lambda_2\)。在 \(t\ge0.1\)，
   \[
    W''+\alpha W'=\sum_{n\ge4}w_n g_n(g_n-\alpha)e^{-g_nt},
      \qquad g_n=\lambda_n-\lambda_2.
   \]
   加权后的主项严格为正，\(t=0.1\) 的尾项/主项上界小于 0.004，
   之后逐项下降。严格裕量和邻域统一谱界保证小参数扰动后仍成立，
   不依赖 \(w_2,w_3\) 的符号或大小。故 \(e^{\alpha t}W'\) 严格递增，
   从而 \(W\) 在该时间区间最多有两个不同零点。
   \([0.002,0.1]\) 的严格正号也由紧性保持。因此此前构造的晚时两零点
   在 \(t\ge0.002\) 恰好是全部零点，不会和另一个这一窗口内的负段叠加。
3. 真实捕获差在整个小邻域的首逆序一侧对所有 \(t>0\) 为正，
   包括尚未验算 W 的极早段。写 \(p=\mu_{R,1}/\mu_{L,1}\)，
   \(\rho=R_L/R_R\)，\(G=F*\nu\)，\(F=A_L\)，\(\nu\) 是其它正常项
   组成的正等待律。精确逆序因子给
   \(
      A_R=\rho[pG-(p-1)(E_{\mu_{R,1}}*G)]\le \rho p F
   \)。
   临界点 \(p=1\)、\(\rho<1\)，故邻域保持 \(\rho p<1\)。这并不证明
   邻域里任意早时的辅助函数 \(W\) 为正；只是这些辅助负段即使存在，
   也不对应物理交叉。

## 方法文献

- Samuel Karlin, *Total positivity, absorption probabilities and applications*
  (1964), [原文](https://www.ams.org/tran/1964-111-01/S0002-9947-1964-0168010-2/S0002-9947-1964-0168010-2.pdf)。
  §0 给出全正核的变号缩减，§5 讨论扩散首次通过。要用在这里，需要对实际
  差值重新构造带有限次变号的输入，并证明具体的时移卷积核满足所需阶数；
  非负性本身不够。不要把对数谱比的核直接当作实际卷积输入。
- I. J. Schoenberg, *On Pólya frequency functions II*
  (1950), [原文](https://acta.bibl.u-szeged.hu/13626/)。
  对时移卷积应采用这类结论，而不是仅说热核在空间变量上全正。
- D. V. Tokarev, *A generalisation of Descartes' rule of signs*,
  [DOI](https://doi.org/10.1017/S144678871100156X)。适用于指数核
  \(e^{-st}\) 的符号测度分块。直接使用本模型全部谱项不够：
  至少在基准中高模权重反复交替，不能由“一次半谱逆序”推出全谱系数只两次变号。
- S. Angenent, *The zero set of a solution of a parabolic equation*
  (1988), [原文](https://math.jhu.edu/~js/Math745/angenent.par1.pdf)。
  控制的是各时间切片的空间零点计数；\(\Delta(t)\) 是固定初始位置下两壁的
  累计差，不能直接把空间零点结论解释成时间上的反超上界。

复算：python first_inversion_weighted_critical_positive.py。
