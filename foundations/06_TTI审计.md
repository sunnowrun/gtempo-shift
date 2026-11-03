# TTI 与均匀局部张量审计

本审计的条件是：固定步长和固定分支 Hamiltonian 下，K、IF 及其本地积分映射均有相同的重复局部张量。相关矩阵 Toeplitz 只是 IF 算法接口的一个条件，不能代替整个网络的检验。

## 1. 结构与实际实现

完整局部指标、费米宇称位、矩阵乘法顺序及折叠公式见推导第 6 节。虚时主体每单元使用 A、B，裸积分后使用同一个 transfer tensor。实时折叠后每层使用同一个作用于矩阵 X 的超张量，并保留反支次数宇称。mixed 使用准备段重复张量与实时重复张量，通过共享浴虚指标的接口连接。

初态、普通迹的反周期闭合、轮廓转折、两个奇插入及淬火是边界/缺陷。纯虚时闭合不等同于无限开链；mixed 不具有物理上不存在的全轮廓单张量平移对称。

| 审计 | 实际比较 | 最大误差 |
|---|---|---|
| combined_uniform_transfer | 重复局部张量 vs 直接算符/全部系数 | 1.24e-14 |
| folded_repeated_K_IF_and_integral | 重复局部张量 vs 直接算符/全部系数 | 3.33e-16 |
| mixed_repeated_segments_quench | 重复局部张量 vs 直接算符/全部系数 | 1.13e-15 |
| mixed_repeated_segments_odd | 重复局部张量 vs 直接算符/全部系数 | 2.69e-15 |

## 2. 原 IF 与上一方案比较

N＝8，δt＝0.1、δτ＝0.1125。残差是相邻行列同时平移一格的最大系数差；实时反支先按物理时间顺序重新排列。

| 方案 | 轮廓 | 同类主体块残差 | 虚实交叉块残差 |
|---|---|---|---|
| order1 | imag | 1.3e-18 | — |
| order1 | real | 1.32e-18 | — |
| order1 | mixed | 1.74e-18 | 0.000195 |
| bath_corrected | imag | 1.74e-18 | — |
| bath_corrected | real | 1.78e-05 | — |
| bath_corrected | mixed | 1.81e-05 | 0.000216 |
| weyl | imag | 1.3e-18 | — |
| weyl | real | 1.32e-18 | — |
| weyl | mixed | 1.74e-18 | 0.000195 |

Weyl 与 order1 的 IF 系数在相同布局下完全相同；实际 Julia 输入数组的逐项相等检查保存在 if_identity.txt。bath_corrected 的实时平移残差不是把容差调小即可消失的截断误差，它来自不同的离散浴目标。mixed 的 cross block 即使原物理浴也不是 Toeplitz；保留这一块及接口是正确结构的一部分。

## 3. uniform/infinite 的界限

已经实现：有限浴、有限链、真正相同的局部张量重复构成的 K/IF 和积分网络，与直接构造逐系数/标量比较。奇插入的 mixed 也已比较，两个奇缺陷的额外符号保留。

未实现：infinite MPS 固定点迭代及收敛判断、无限链环境、连续浴有限键维数定理。实际 package 的有限链 canonicalize 后站点张量不要求数值逐项相等；重复性成立于明确给出的局部费米网络，不应把有限规范化结果误称为 uniform MPS 求解器。

现有 IF 构造兼容范围：PartialIF＋zipup 三轮廓实际执行；TTI-IF/ExactTTI 虚时及实时的原始运行结果见 julia_tti.csv；mixed 专用 TTI 软件不存在，未声称运行。
