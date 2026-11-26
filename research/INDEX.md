# 三层介质 Fick 扩散：累计捕获与反超

单位区间、中点释放、理想介质界面，两壁具有同一 Robin 反应性。主要问题：最终分配与短时到达顺序不能决定有限时间的全部捕获顺序，三层介质究竟能反超几次？

**已确证：** 存在恰好三次反超的参数；共同壁面反应性可使一对穿越产生或消失；存在一个临界点附近，恰好三次反超可全部放在任意给定截止之后。另一个临界点附近，一次早期反超可以和任意晚的一对反超共存。各结果的局部上限不外推为全三层模型族定理。未发现超过三次的已认证例子。

## 证明与证书

1. [模型、早期与最终捕获量、三个严格符号反转的下界](README.md)：有理区间变换验证；另有数值反演作诊断。
2. [固定算例恰好三次](EXACTLY_THREE.md)：全时间区间证明与根的严格夹区。
3. [有限时间双零点与最终打平](BIFURCATIONS.md)：只调共同反应性可产生或消灭一对反超。
4. [一次早期＋两次晚期](LATE_PAIR_COUNT.md)：控制到无限时间的局部上限。
5. [三次全部任意晚](LATE_TRIPLE.md)：两个慢模同时消失、三个局部坐标与加权 Rolle 上界。
6. [临界族解析降维与数值延续](CRITICAL_BRANCH.md)：单点严格证书之外的探索。结果不当作证明。

证书涉及解析短时估计、本征值与谱尾误差界、区间算术及零点计数。普通浮点只负责找候选括区，或单独提供探索数据。这里是 Fick 模型研究，不是已投 JSP 的自然尺度双 Robin 扩散稿。

## 复跑

验证环境：Python 3.12.14，NumPy 2.3.5，SciPy 1.17.0，mpmath 1.3.0。进入本目录：

    python verify.py
    OPENBLAS_NUM_THREADS=1 python upper.py
    OPENBLAS_NUM_THREADS=1 python certify_folds.py
    OPENBLAS_NUM_THREADS=1 python certify_late_pair.py
    OPENBLAS_NUM_THREADS=1 python late_triple.py

相应运行输出见本目录的文本文件。探索脚本 `critical_branch.py` 和 `scan_critical.py` 不参与严格结论。没有为这个归档增加 GitHub Actions。
