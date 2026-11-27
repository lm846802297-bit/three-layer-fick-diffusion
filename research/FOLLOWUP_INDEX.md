# 三层 Fick 扩散：后续研究归档

[基础模型与复跑说明](INDEX.md)

## 后续归档：首谱逆序、零点与方法限制

- [慢模消失的普遍必要性和局部构造](LATE_CROSSING_THEOREM.md)
- [半谱排序与可严格排除的参数区](GLOBAL_SIGN_REDUCTION.md)；[带符号延迟](SIGNED_LOW_MODES.md)；[单逆序降维](ONE_INVERSION_REDUCTION.md)；[凸谱比较](CONVEX_DELAY_BRIDGE.md)
- [恰有一个逆序时恰好三次反超的实例](ONE_MODE_THREE_REVERSALS.md)
- [首模逆序的卷积障碍](FIRST_INVERSION_OBSTRUCTIONS.md)；[任意晚的实际反超对](FIRST_INVERSION_LATE_PAIR.md)；[有限时间双零点](FIRST_INVERSION_FINITE_FOLD.md)
- [全时间正性判据](FIRST_INVERSION_RESIDUAL_BAND.md)；[负谷的凸包络](FIRST_INVERSION_VALLEY_ENVELOPE.md)
- [带权零点局部计数](FIRST_INVERSION_WEIGHTED_COUNT.md)；[改变几何后的复核](FIRST_INVERSION_WEIGHTED_GEOMETRY.md)；[带权函数相切检验](FIRST_INVERSION_WEIGHTED_TANGENCY.md)
- [带权导数两个任意晚零点](FIRST_INVERSION_WEIGHTED_LATE_PAIR.md)；[双消失点的区间正性](FIRST_INVERSION_WEIGHTED_CRITICAL_POSITIVE.md)
- [全正核路线与谱符号法的局限](FIRST_INVERSION_VARIATION_DIMINISHING.md)
- [物理输入 J 的任意晚双零点](FIRST_INVERSION_PHYSICAL_INPUT_LATE_PAIR.md)：辅助函数变号不等于实际反超。
- [晚时 J 零点对不能叠加早期过零](FIRST_INVERSION_PHYSICAL_INPUT_NO_STACKING.md)：规则、非均匀临界极限的解析结果。
- [单个晚时 J 零点与静态打平面](FIRST_INVERSION_SINGLE_LATE_ZERO.md)：严格逆序时的唯一晚时上穿率，以及双消失邻域真实捕获差的全时间正性。
- [静态打平面的实 Laplace 变换](FIRST_INVERSION_STATIC_TIE_LAPLACE.md)：变换恰有一个简单零点；显式反例说明这一性质不能反推出时间域零点上限。
- [静态打平面的时间域单零点判据](FIRST_INVERSION_STATIC_TIE_TIME.md)：消去逆序极点后的公共核严格对数凹，给出仅用正常半谱的风险率判据；临界邻域严格只有一个简单零点。
- [均匀右半区的风险率晚时界](FIRST_INVERSION_STATIC_TIE_UNIFORM_RIGHT.md)：Robin 特征函数精确消去高频积，给风险率阈值以上的显式截止；中间窗仍未关闭。
- [固定打平基准的全时间单零点证书](FIRST_INVERSION_STATIC_TIE_MIDDLE_CERTIFICATE.md)：严格认证中间窗和两个端点段，给该基准唯一风险率上穿与唯一简单 J 零点。
- [静态打平面的剥离输入单零点](FIRST_INVERSION_STATIC_TIE_STRIPPED_INPUT.md)：方差张量化与右谱级数控制去掉共同高阶等待后的单下穿；未证明传到完整 J。
- [静态打平首逆序的全时间单零点定理](FIRST_INVERSION_STATIC_TIE_COMPLETE.md)：逐段恢复高阶指数等待，封闭打平子域及有严格光学间隙的左壁最终落后侧；不覆盖整个模型。
- [同一几何的真实反超对照](FIRST_INVERSION_STATIC_TIE_GROUND_CONTRAST.md)：κ=4 与 8 的 J 都只有一个零点，但累计差分别反超一次和零次。
- [固定几何的 Robin 反超临界值](FIRST_INVERSION_STATIC_TIE_ROBIN_TRANSITION.md)：整个 κ∈[1,10] 的时域复核与实际反超从一次变零次的唯一临界 κ。
- [左壁最终占优侧的局部折叠幂律](FIRST_INVERSION_POSITIVE_STATIC_FOLD_LAW.md)：分类静态打平附近的真实反超，并证明已认证临界几何有唯一折叠曲线。

附录脚本和原始运行输出在本目录。部分脚本为浮点筛查，不能作为全域证明。后续归档完整性可用 SHA-256 校验：
    
    sha256sum -c FOLLOWUP_SHA256SUMS

**首模逆序子域以及整个三层模型的全时间反超上限均仍未证明。**
