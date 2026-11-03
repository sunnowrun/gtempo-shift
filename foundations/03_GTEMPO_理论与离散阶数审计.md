# GTEMPO 的离散目标、张量网络表示与阶数审计

> **2026-10-21 更正：本文涉及源重建的一般推导与一般一阶定理尚未完成核验。** 撤回把第 4.1—4.3 节及相关概括当作已建立严格证明的说法。普通相干态恒等分解与候选源重建分解必须区分；固定源的浴 Gaussian 恒等式与实际 QuAPI 核离散化及共同积分之间的等价性，还需完整论证。有限模型脚本保留为指定离散公式的计算记录，不能据此替代一般证明或声称已直接复现 Julia package。此更正不主张上述候选重建必然错误，也不主张 GTEMPO 不可推导。

日期：2026-10-20。本文只在当前工作目录创建文件；论文与 package 均按只读参考处理。

## 1. 可以下什么结论

**可以建立坚实推导，但必须区分三种命题。**

1. 给定一个有限时间格点上的 Grassmann 积分，将其表示为 GMPS 并收缩，是可以严格证明的代数事实。未截断时，这一步不带时间离散误差，也不决定算法是一阶还是二阶。
2. 从真实 Hamiltonian 的 Lie 或 Strang 离散传播子出发，插入相干态恒等式、再严格积掉 Gaussian 浴，可以严格得到某个离散 $K\mathcal I$。这条路线具有清楚的一阶或二阶误差定理。
3. **把连续 IF 的核按小矩形积分，再与当前 $K$ 相乘，不自动等于上一项得到的离散 IF。** 原始论文中将 Trotter 阶数直接赋给最终 GTEMPO 计算，需要补上这个等价性检查。

本文对当前源码对应的一个可解子类给出完整反例：杂质自由 Hamiltonian 为零、只有一个零能浴模时，所谓 `order=2` 的未截断离散目标仍有一般非零的线性步长误差。此时论文写出的 Strang 分解本身完全精确，因此不能把这个误差归因于 Strang 的交换子。

这说明存在具体的离散化缺口；它不意味着 Grassmann 张量网络方法没有理论基础。

## 2. 对照材料与证据边界

主要核查：

- [原始 GTEMPO 论文检索副本 GT013](/home/zjsun/Desktop/Agent/my_paper/thesis2/gtempo-refs/markdown/GT013.md)：Phys. Rev. B **109**, 045140，DOI 10.1103/PhysRevB.109.045140；正文式 (3)、(4) 给出 Lie/Strang 分解，附录 A 式 (A4) 给出核的矩形积分。这里区分论文声明与下文独立推导。
- [GTEMPO 综述 GT023](/home/zjsun/Desktop/Agent/my_paper/thesis2/gtempo-refs/markdown/GT023.md)：J. Chem. Phys.，DOI 10.1063/5.0226167；式 (43)—(48) 分别涉及核离散化、$K$ 与系数张量的 MPS 表示。
- [correlationfunction.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/correlationfunction.jl:11)：`order=2` 调用双倍时间格点的浴核。
- [realtime_stepper.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/influencefunctional/hybridization/partialif/realtime_stepper.jl:81)：`top_right_if_final` 与 `top_right_if` 合并半步核；最终观测必须使用 `finalize=true` 的 IF。
- [QuAPI/fermionic.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/realtime/fermionic.jl:92) 与 [QuAPI/util.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/realtime/util.jl:57)：给出实际同格点、异格点系数。
- [QuAPI/corrmat.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/corrmat.jl:19)：对角元是两个存储向量首项之和，不能只读其中一项。
- [siam.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/predefined/siam.jl:237)：系统步进核；本例取参数为零，仅留下相干态重叠。
- [boundarycondition.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/gvconnections/boundarycondition.jl:54) 与 [gf.jl](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/gf.jl:16)：迹边界及归一化 Grassmann 关联函数。

数值验证是 **Python 对以上源码公式的独立重建，没有执行 Julia package**。它验证离散数学目标，而不是 package 的所有运行路径。单浴模的谱为 delta 谱，不依赖当前谱函数 API 能否直接输入 delta 谱；也可直接构造有限模的关联矩阵。本例的浴核积分解析完成，没有谱积分误差。

原始 PDF 的提取文本有公式和栏序损坏；本文不以损坏公式作为证明，所有反例公式均从源码及有限维模型重建。源文件 SHA-256 见 [source_manifest.json](source_manifest.json)。

## 3. 首先定义一个真正可控的离散目标

### 3.1 对象、假设、误差参数

考虑有限个杂质轨道和有限个浴模：

$$
H=A+B,\qquad A=H_{\rm imp},\qquad
B=H_{\rm bath}+H_{\rm hyb},
$$

$$
H_{\rm bath}=\sum_k\epsilon_k b_k^\dagger b_k,\qquad
H_{\rm hyb}=\sum_{\alpha k}
\bigl(V_{\alpha k}d_\alpha^\dagger b_k+
V_{\alpha k}^* b_k^\dagger d_\alpha\bigr).
$$

初态取杂质与 Gaussian 浴的乘积。杂质可有局域相互作用。固定物理时间 $t=Nh$，研究 $h\to0$，并先排除 MPS 截断、记忆截断及谱积分误差。

有限浴假设使以下算符范数论证无歧义；把它推广到连续、无界谱浴，需要另外控制谱矩和相应的交换子或 Duhamel 余项，不能把有限维常数自动当成浴模数无关的常数。

### 3.2 一阶来自局域算符误差与稳定性

定义

$$
U_{1,h}=e^{-ihB}e^{-ihA}.
$$

直接展开可得

$$
U_{1,h}-e^{-ih(A+B)}
=\frac{h^2}{2}[A,B]+O(h^3).
$$

利用两种传播子的幺正性，逐项望远镜相减：

$$
\left\|U_{1,h}^N-e^{-itH}\right\|
\leq N\left\|U_{1,h}-e^{-ihH}\right\|
\leq C_t h.
$$

归一化初态的约化态与有界观测量也继承这一上界。例如

$$
|\langle O\rangle_h-\langle O\rangle|
\leq 2\|O\|\left\|U_{1,h}^N-e^{-itH}\right\|.
$$

这里 $O(h)$ 是普遍上界，不表示每个初态、每个观测量必有非零线性系数。对易模型和对称性可能使某些误差消失。

Strang 步

$$
U_{2,h}=e^{-ihB/2}e^{-ihA}e^{-ihB/2}
$$

满足局域误差 $O(h^3)$，因此固定时间的全局误差是 $O(h^2)$。

**这些定理能传给张量网络计算的前提，是网络精确表示相应的离散传播子。网络名字或时间格点布局不能代替这个前提。**

### 3.3 从离散传播子到相干态积分

固定相干态约定与 Grassmann 积分次序，使

$$
1=\int d\bar a\,da\,e^{-\bar a a}|a\rangle\langle\bar a|,
\qquad \langle\bar a'|a\rangle=e^{\bar a'a}.
$$

在每个离散因子之间插入恒等式，保留杂质矩阵元

$$
K_h=\prod_{\text{contour links}}
\langle\bar a_{\rm out}|e^{\mp ihA}|a_{\rm in}\rangle
\times(\text{initial and trace factors}).
$$

对相互作用杂质，不应未经检查就把这个精确矩阵元换成 Hamiltonian 正常符号的简单指数。

因 $B$ 是二次型，设其单粒子矩阵为 $\mathsf B$，则精确一步的相干态核是

$$
\langle\bar z'|e^{-ihB}|z\rangle
=\exp\bigl(\bar z'e^{-ih\mathsf B}z\bigr),\qquad z=(a,b).
$$

将 $e^{-ih\mathsf B}$ 分为杂质—浴四个块，可以逐项读出杂质局域修正、杂质—浴交叉项和浴传播块。Strang 端点用半步矩阵。

把所有浴变量收集后，积分具有以下 Gaussian 结构：

$$
\int D[\bar b,b]\,
\exp\bigl[-\bar b M_Bb+\bar a Cb+\bar b Da\bigr]
=\det(M_B)\exp\bigl(\bar a CM_B^{-1}Da\bigr).
$$

还应计入上一步单粒子矩阵中直接作用于杂质的局域块。归一化浴迹后，得到

$$
\mathcal I_{p,h}^{\rm exact\ discrete}
=c_h\exp\bigl(\bar a\,\Gamma_{p,h}\,a\bigr),
\qquad p=1,2.
$$

因此

$$
\langle O\rangle_{p,h}
=\frac{\int D[\bar a,a]\,
\mathcal O_{p,h}K_{p,h}\mathcal I_{p,h}^{\rm exact\ discrete}}
{\int D[\bar a,a]\,K_{p,h}\mathcal I_{p,h}^{\rm exact\ discrete}}
$$

是有明确误差阶数的离散计算目标。算符插入也必须对应同一个离散轮廓。

这一推导对杂质相互作用的大小没有弱耦合假设。它也不要求浴 IF 的张量网络一定容易压缩。

## 4. 当前核积分步骤为何需要额外证明

连续时间浴积分给出形式上的精确 IF：

$$
\mathcal I[\bar a,a]
=\exp\left[-\int_C ds\int_C ds'\,
\bar a(s)\Delta(s,s')a(s')\right].
$$

当前常见做法选择分片常数源，并定义

$$
\eta_{ij}^{\zeta\zeta'}
=-\int_{C_i^\zeta}ds\int_{C_j^{\zeta'}}ds'\,
\Delta^{\zeta\zeta'}(s,s'),\qquad
\mathcal I_h=\exp\left(\sum_{ij\zeta\zeta'}
\bar a_i^\zeta\eta_{ij}^{\zeta\zeta'}a_j^{\zeta'}\right).
$$

这里把负号吸收到 $\eta$ 中，与 QuAPI 源码约定一致。

对于固定的外部 Grassmann 源，这种矩形积分有明确意义，也可以精确计算。但是最终还要将这些源用相干态测度和 $K_h$ 积掉。**源的格点识别、相干态重叠以及同格点接触项，都属于离散化的一部分。**

因此以下两条路线不能直接交换：

$$
\text{精确离散传播子}\longrightarrow\text{相干态积分}
\longrightarrow\text{积掉浴}
$$

和

$$
\text{连续 IF}\longrightarrow\text{分片常数源与核积分}
\longrightarrow\text{与另行构造的 }K_h\text{ 相乘并积分}.
$$

若要证明二阶，必须证明两者的**收缩结果**之差至多为 $O(h^2)$，或证明一个等价稳定步进器的局域误差为 $O(h^3)$。不能只对核的标量求积精度做论证，也不能假设 Grassmann 路径是一条可对时间做普通 Taylor 展开的光滑数值函数。

下面的反例直接在最终积分层面检查这一点。

### 4.1 将分片常数源重建为局域步进器

对于有限浴，还能补上一般的一阶上界。这里采用按完整矩形定义、并保留同分支时间排序的 IF；这与下面核查的 QuAPI 系数一致。

令 $F_\alpha=\sum_kV_{\alpha k}b_k$，在一个宽度为 $w$ 的单元中，把杂质变量作为外部 Grassmann 源，定义浴算符

$$
\mathscr B(\bar a,a)=H_{\rm bath}
+\sum_\alpha(\bar a_\alpha F_\alpha+F_\alpha^\dagger a_\alpha),
\qquad E_w(\bar a,a)=e^{-iw\mathscr B(\bar a,a)}.
$$

所有外部 Grassmann 奇元与浴费米奇算符反对易。浴在单元内的时间仍连续，源在单元内保持不变。逐单元相乘 $E_w$ 再取 Gaussian 浴迹，用 Wick 定理积掉浴，恰好给出这些分片源的矩形积分 IF。这一步不是近似截断杂化展开。

接着对一个杂质节点使用相干态分解。定义带正确分次次序的线性映射

$$
\mathscr Q[f]=\int D[\bar a,a]\,
e^{-\bar a a}|a\rangle f(\bar a,a)\langle\bar a|,
\qquad R_w=\mathscr Q[E_w].
$$

这是把对角相干态源重新变成算符；需要保留浴奇算符与源的相对次序。在一致相干态约定下，单模恒等式为

$$
\mathscr Q[1]=1,\qquad
\mathscr Q[a]=d,\qquad
\mathscr Q[\bar a]=d^\dagger,
$$

$$
\mathscr Q[a\bar a]=dd^\dagger=1-d^\dagger d,\qquad
\mathscr Q[\bar a a]=d^\dagger d-1.
$$

前三个等号可用相干态的本征关系和恒等分解逐项验证；后两个等号还用到 Grassmann 反对易关系。由此可见 $\mathscr Q$ 是线性映射，**不是保乘法的映射**。特别是

$$
\mathscr Q[\bar a a]\neq
\mathscr Q[\bar a]\mathscr Q[a].
$$

将完整相干态路径积分按节点积分，前向传播子正是

$$
P_h=R_{w_N}e^{-ihA}R_{w_{N-1}}e^{-ihA}\cdots
e^{-ihA}R_{w_0}.
$$

该等价关系也可反过来证明：把每个 $R_w$ 的积分定义代入右边，杂质核和分片源浴算符逐项恢复；再取浴迹就恢复 $K_h\mathcal I_h$。

短时展开给出

$$
R_w=1-iwB-\frac{w^2}{2}\mathscr Q[\mathscr B^2]+O(w^3),
$$

$$
R_w-e^{-iwB}
=-\frac{w^2}{2}\left(\mathscr Q[\mathscr B^2]-B^2\right)+O(w^3).
$$

一阶项正确，是因为耦合对源线性；二阶项一般不正确，是因为上述映射不保乘法。这明确指出了“连续 IF 分片后再与 $K$ 积分”和“先离散真实传播子再积掉浴”的差别。

### 4.2 有限浴的一阶误差上界

固定有限浴、有限杂质轨道数、有界 $A$ 和固定 $t$。上面的所有 Grassmann 空间及算符空间都有限维，因此存在与 $h$ 无关的局域常数 $C$，使

$$
\|R_w-e^{-iwB}\|\leq Cw^2,\qquad
\|R_w\|\leq1+C|w|
$$

对足够小的 $w$ 成立。半步端点布局满足

$$
\sum_j w_j=t,\qquad
\sum_j w_j^2=th-h^2/2.
$$

在 $P_h$ 中逐个用真实浴—耦合传播子替换 $R_{w_j}$，乘积范数由 $e^{Ct}$ 控制，故

$$
\|P_h-U_{2,h}^N\|
\leq C_te^{Ct}\sum_jw_j^2=O(h).
$$

再加上真实 Strang 传播子的 $O(h^2)$ 误差，得到

$$
\|P_h-e^{-itH}\|=O(h).
$$

前后向传播、浴迹与固定数目的有界算符插入均保持这个阶数。由于 $R_w$ 不保证幺正，离散分母未必严格为 1，但它满足 $Z_h=1+O(h)$，足够小步长时有界远离零；除以 $Z_h$ 仍保持 $O(h)$ 上界。整宽端点布局则额外有 $O(h)$ 的端点传播，结论相同。

这构成了**有限浴、完整矩形分片源 IF、精确系统核、无压缩误差时的一般一阶上界**。下一节的反例证明其中的线性误差确实可以非零，因而不能普遍提升为二阶。

当前 QuAPI 的跨分支对角系数也须按实际矩阵索引规则核查，见第 8 节；不能把内部存储的半个三角形误当成最终完整矩阵元。

### 4.3 单浴模中可以直接看见接触项

对下一节的 $H_{\rm bath}=0$ 单模模型，

$$
\mathscr B=v(\bar a b+b^\dagger a),\qquad
\mathscr B^2=v^2\bar a a(1-2n_b),\qquad
\mathscr B^3=0.
$$

因此上述源重建的局域算符精确为

$$
R_w=1-iwB-\frac{v^2w^2}{2}(n_d-1)(1-2n_b).
$$

真实 $B^2=v^2(n_d+n_b-2n_dn_b)$，于是

$$
R_w-e^{-iwB}=\frac{v^2w^2}{2}(1-n_b)+O(w^3).
$$

这里的额外项就是不能把相干态符号相乘当成算符相乘的直接后果。真空矩阵元为 $1+v^2w^2/2$。在单粒子子空间 $\{|d\rangle,|b\rangle\}$ 中除以这个真空矩阵元，再把浴幅度写成 $-iy$，正好得到下一节从 Gaussian 积分导出的实矩阵 $M(w)$。两种重建相互吻合。

## 5. 一个没有任何 Trotter 误差的反例

### 5.1 模型与严格答案

取一个无自旋杂质模 $d$ 和一个浴模 $b$：

$$
H_{\rm imp}=H_{\rm bath}=0,\qquad
H=v(d^\dagger b+b^\dagger d),\qquad v\in\mathbb R.
$$

初态为两个模的真空，计算

$$
g(t)=\langle d(t)d^\dagger(0)\rangle=\cos(vt).
$$

如果采用通常带 $-i$ 的 greater Green 函数定义，只需整体乘以 $-i$，误差阶数不变。

因为 $A=0$，Lie 和 Strang 的真实算符乘积都严格等于 $e^{-itH}$，对任意步长成立。这里没有 Hamiltonian 分裂误差，没有浴截断，没有记忆截断。

### 5.2 当前 IF 与 K 给出的 Gaussian 矩阵

节点记作 $j=0,\ldots,N$。对源码标为二阶、并正确 `finalize=true` 的目标，半步核合并给出有效宽度

$$
w_0=w_N=h/2,\qquad w_1=\cdots=w_{N-1}=h,\qquad Nh=t.
$$

前向分支的核为

$$
\eta^{++}_{ij}=
\begin{cases}
-v^2w_iw_j,&i>j,\\
-v^2w_i^2/2,&i=j,\\
0,&i<j.
\end{cases}
$$

同格点的 $1/2$ 来自前向分支时间排序后的三角形积分，而不是人为近似掉浴频率。本例的核在时间上为常数，这个三角形积分本身是精确的。

系统核只有

$$
K_+=\exp\left(\sum_{i=1}^N\bar a_i a_{i-1}\right).
$$

将测度 $e^{-\sum_i\bar a_i a_i}$ 一并计入，前向 Gaussian 矩阵 $Q$ 满足

$$
Q_{ii}=1+v^2w_i^2/2,\qquad
Q_{ij}=-\delta_{i,j+1}+v^2w_iw_j\quad(i>j),\qquad
Q_{ij}=0\quad(i<j).
$$

Gaussian 恒等式给出归一化关联函数

$$
g_h(t)=(Q^{-1})_{N0}.
$$

这并非丢掉了后向分支的假设：在本例真空初态下 $\eta^{+-}=0$，初态相干态矩阵元为 1。把前向节点、迹变量、后向节点依次排列，完整 Keldysh 矩阵为块下三角矩阵。前向—前向逆矩阵块正是 $Q^{-1}$。脚本另外重建了完整矩阵验证这一结论。

### 5.3 精确消元成 2×2 递推

令 $x_i=(Q^{-1})_{i0}$，则

$$
(1+v^2w_0^2/2)x_0=1,
$$

$$
(1+v^2w_i^2/2)x_i
=x_{i-1}-v^2w_i\sum_{j<i}w_jx_j\quad(i\geq1).
$$

定义 $y_i=v\sum_{j=0}^i w_jx_j$，并取辅助初值 $x_{-1}=1,y_{-1}=0$。所有节点统一写成

$$
\begin{pmatrix}x_i\\y_i\end{pmatrix}
=M(w_i)\begin{pmatrix}x_{i-1}\\y_{i-1}\end{pmatrix},
$$

$$
M(w)=\frac{1}{1+v^2w^2/2}
\begin{pmatrix}
1&-vw\\
vw&1-v^2w^2/2
\end{pmatrix}.
$$

这是对当前离散 Grassmann 积分的**精确恒等变换**，没有再次截断时间展开。

于是所谓二阶目标恰好等于

$$
g_h(t)=e_1^T M(h/2)M(h)^{N-1}M(h/2)e_1,
\qquad e_1=(1,0)^T.
$$

这也是一个特别容易独立重现的闭式离散目标。

### 5.4 局域缺陷与固定时间的线性误差

令

$$
L=\begin{pmatrix}0&-v\\v&0\end{pmatrix},\qquad
D=\begin{pmatrix}0&0\\0&-v^2/2\end{pmatrix}.
$$

展开精确矩阵 $M(w)$ 得

$$
M(w)=I+wL+w^2\begin{pmatrix}-v^2/2&0\\0&-v^2\end{pmatrix}+O(w^3),
$$

$$
\log M(w)=wL+w^2D+O(w^3).
$$

真实传播的辅助二维表示为 $e^{wL}$。因此

$$
M(w)-e^{wL}=w^2D+O(w^3).
$$

**这就是一阶的局域来源：不是核积分算错，而是最终积分对应的步进器带有 $O(w^2)$ 的额外项。**

左右两个半步不会抵消内部 $N-1$ 个这样的缺陷。利用有限维矩阵指数的一阶扰动展开，固定 $t$ 时

$$
g_h(t)=\cos(vt)+h\,e_1^T\left[
\int_0^t e^{(t-s)L}D e^{sL}\,ds\right]e_1+O(h^2).
$$

这里可以显式算出积分：

$$
e_1^T e^{(t-s)L}D e^{sL}e_1
=\frac{v^2}{2}\sin(v(t-s))\sin(vs),
$$

所以

$$
\boxed{\displaystyle
g_h(t)-\cos(vt)
=\frac{h}{4}\left[v\sin(vt)-v^2t\cos(vt)\right]+O(h^2).}
$$

括号内一般非零，故这是严格的一阶反例，而不仅是误差曲线看起来像一阶。在系数偶然为零的特定时间点，单独测该点可能看不出一阶，不能据此推断整个方法二阶。

还可检查

$$
M(-h)M(h)=I+2h^2D+O(h^3).
$$

因此当前有效步进器本身不具备对称传播公式所需的逆时间性质。“节点左右对称”不足以让它成为真正的 Strang 方法。

### 5.5 一阶标签的端点还产生另一项误差

当前 `order=1` 使用 $N+1$ 个整宽节点，故在这个模型中

$$
g_{1,h}(t)=e_1^TM(h)^{N+1}e_1.
$$

其有效源区间总宽度为 $t+h$，于是

$$
g_{1,h}(t)-\cos(vt)
=h\left[-v\sin(vt)
+\frac{v\sin(vt)-v^2t\cos(vt)}4\right]+O(h^2).
$$

二阶标签的半步端点处理消除了这里的第一项，却留下上一节的体内线性项。这解释了为什么它可以明显更准确，同时仍然只有一阶。

## 6. 独立数值检查

运行：

```bash
python verify_discretization.py
```

参数 $v=0.7,t=1$，严格答案为 $0.764842187284$。表中误差均针对这个固定物理时间。

| 源码标签 | N | 绝对误差 | 相邻加密的测量阶数 |
|---|---:|---:|---:|
| order=1 | 40 | 0.01086758449 | 1.00874 |
| order=1 | 80 | 0.005416557158 | 1.00458 |
| order=1 | 160 | 0.002703879946 | 1.00235 |
| order=1 | 320 | 0.001350829111 | 1.00119 |
| order=2 | 40 | 0.0004865552804 | 1.02857 |
| order=2 | 80 | 0.0002407131267 | 1.01529 |
| order=2 | 160 | 0.0001196991052 | 1.00790 |
| order=2 | 320 | 0.00005968315262 | 1.00402 |

阶数定义为

$$
p_h=\log_2\frac{|g_{2h}(t)-g(t)|}{|g_h(t)-g(t)|}.
$$

解析预测的线性系数分别为

$$
c_1=-0.431907453742,\qquad c_2=0.0190449273242.
$$

脚本验证了误差除以步长趋于这两个数。它还互相核对完整 Keldysh Gaussian 积分、前向三角求解、2×2 递推和源重建的 4×4 Fock 空间算符乘积，最大差异为 $8.88\times10^{-16}$；另有一个小规模 Grassmann 多项式直接 Berezin 积分检查。完整数据见 [convergence.csv](convergence.csv)。

## 7. 为什么张量网络表示本身是严格的

固定全局变量顺序 $\xi_1,\ldots,\xi_L$，任意有限 Grassmann 多项式唯一写为

$$
F(\xi)=\sum_{n_1,\ldots,n_L\in\{0,1\}}
F_{n_1\cdots n_L}\xi_1^{n_1}\cdots\xi_L^{n_L}.
$$

系数就是普通有限维张量。对连续切口做精确 SVD，可得到精确 MPS：

$$
F_{n_1\cdots n_L}
=A_1^{n_1}\cdots A_L^{n_L}.
$$

这个存在性不保证小键维度；一般最坏情形需要指数大的键维度。

Grassmann 乘法对应系数卷积和明确的交换符号。若两项的占据向量是 $p,q$，有公共变量时乘积为零；无公共变量时重排的符号为

$$
(-1)^{\sum_{i>j}p_iq_j}.
$$

此外 IF 的每个二次项是偶元，彼此对易，因此

$$
\mathcal I_h
=\prod_i\exp\left(\bar a_i\sum_j\eta_{ij}a_j\right)
=\prod_i\left(1+\bar a_i\sum_j\eta_{ij}a_j\right).
$$

最后一个等号来自行因子的平方为零。这正是 partial-IF 构造的严格代数依据，不是另一次一阶 Trotter 分解。GMPS 相乘中加入的费米符号负责实现上面的卷积符号。

最后，Berezin 积分是系数空间上的线性泛函；相干态测度也是已知的有限多项式。把它们逐站点收缩就是普通张量网络的线性收缩。

因此，“IF 与 K 构成 Grassmann 多项式 → GMPS → 积分”的链条可以完全严格。截断键维度会额外引入压缩误差，但不会修复错误的时间离散目标。

## 8. 跨分支对角元的源码一致性检查

仅阅读 `fermionic.jl` 第 118、132 行，会看到跨分支同物理格点的存储项也使用三角积分；但这不是最后的矩阵元。`CorrelationMatrix.getindex` 在对角位置将两个存储向量的首项相加，见 `corrmat.jl` 第 19 行。

例如本例空浴的连续跨分支核是常数。单个宽度 $w$ 的节点应按完整矩形积分给出

$$
\eta^{-+}_{ii}=-v^2w^2,
$$

源码恰好把两个 $-v^2w^2/2$ 相加，恢复 $-v^2w^2$。二阶的内节点合并两段半步后，该跨分支系数为 $-v^2h^2$，同样正确。

同分支的对角矩阵元则是占据与未占据两种时间排序三角积分的组合；在本例空浴下，只有其中一个有贡献，保留 $1/2$。验证脚本据此重建了完整四块轮廓核。

因此，本例的线性误差不能归因于漏加跨分支对角系数。它来自第 4、5 节指出的源重建与真实传播之间的二阶局域缺陷。

## 9. 什么已经证明，什么还不能声称

已证明：

- 有限 Grassmann 离散目标的精确张量网络表示与积分有严格基础。
- 对有限浴，从真实 Lie/Strang 步出发严格积掉浴，可以构造有一阶/二阶误差上界的离散 $K\mathcal I$。
- 对有限浴的分片常数源目标，可通过相干态源重建的局域映射证明一般 $O(h)$ 上界；这里显式展示了局域映射和稳定性估计。
- 对当前源码公式对应的上述有限浴模型，`order=1` 与 `order=2` 目标均一般为一阶；后者的线性系数已解析算出。
- 原始论文展示的 Hamiltonian 对称分解阶数，不足以证明当前核积分与 $K$ 组合的二阶精度。

尚未证明或尚未执行：

- 没有证明当前所有 GTEMPO 变体、所有初态、所有连续浴都是一阶；一个反例否定的是普遍二阶声明。
- 没有给出任意无界连续浴下、对谱截止均匀的误差定理。
- 没有在 Julia 中直接运行这个 package，也没有复现用户此前测试的具体参数。
- 没有审计无限时间、虚时间、BCS 浴、所有电流算符插入或特定压缩方法的阶数。

第 4 节的一般上界依赖有限浴或其他能保证这些局域映射与余项一致有界的条件。它不是对任意无界连续浴的统一定理，也不能保证所有特殊观测量恰好一阶。用于博士论文时，应明确列出这些条件，再独立处理希望采用的谱函数及连续浴极限。

## 10. 如何建立真正可审计的二阶版本

最直接的路径是第 3 节：先固定真实 Strang 传播子，精确计算 Gaussian 的 $B$ 半步和整步单粒子矩阵，再通过浴 Gaussian 消元得到离散 IF，最后用同一套 GMPS 代数来表示它。

这条路线中，二阶来自真实传播子的局域误差 $O(h^3)$；浴积分与张量网络代数未截断时都是恒等变换。对连续浴，需要另证明相应极限或导出无需显式有限浴近似的离散核。

最基本的判据是本报告的 $H_{\rm imp}=0$ 模型：若 IF 严格对应论文所写 $B$ 步，这个模型在任意步长都应精确。当前公式没有通过它。先通过这个判据，再在非对易模型中测得固定时间的二阶收敛，才能把“二阶 GTEMPO”作为有理论与数值共同支持的结论。

对论文写作，现阶段比较稳妥的表述是：**Grassmann 张量网络收缩对给定离散作用量是精确表示；常用核积分离散化的时间阶数需要独立分析。当前所核查的对称格点方案存在一阶反例，不能仅按 Trotter 分解名称声称二阶。**
