# 当前 package 中 order=1 实时 GTEMPO：从 Hamiltonian 到代码离散目标

日期：2026-10-21。本文重新核对代码和代数，只讨论 `RealGrassmannLattice1Order`。原有文档不作为证明前提，也不作修改。

本文的主线是预定义 `AndersonIM`（SIAM）的精确系统核、有限正规自由费米浴、因子化初态、杂质末时刻偶观测量。辅助重建映射用于分析代码公式，并不是 package 实际执行的算符步骤。全文区分三件事：有误差控制地逼近真实演化；精确表示这个近似算符；核对它与代码构造的多项式一致。

所有“代码等于”均指精确系数、保留完整记忆、不作 MPS 或积分截断的离散计算目标。浮点、谱求积、初态系数过滤和张量压缩误差另计。本文先在固定有限浴中闭合证明，不声称已经证明任意连续谱或无限带宽极限。

## 1. 首先确定代码实际构造什么

### 1.1 时间索引：系统间隔和源节点不是同一个计数

设代码参数 `δt=h`、`N=N`、`bands=m`。`lattice.k=N+1`、`lattice.t=Nh`、`lattice.Nt=N`。这些定义位于 [realtime.jl:47](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/lattices/realtime.jl:47)。

本文的节点编号从零开始，代码的源节点编号从一开始。对应关系如下。

| 对象 | 本文编号或数目 | 代码编号或位置 |
|---|---|---|
| 系统传播间隔 | $j=0,\ldots,N-1$，共 $N$ 个 | 系统核循环 `i=1:N` |
| 源节点 | $j=0,\ldots,N$，共 $N+1$ 个 | `index(lattice,j+1,...)` |
| 系统标记时间 | $t_j=jh$ | `lattice.ts` |
| 末节点 | $j=N$ | `lattice.k=N+1` |
| 辅助迹变量 | $\chi_\alpha,\bar\chi_\alpha$，共 $m$ 对 | `index(lattice,0,...)`，不属于任何分支 |

因此共有 $4m(N+1)+2m$ 个 Grassmann 生成元，与 [realtime.jl:64](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/lattices/realtime.jl:64) 的长度公式一致。代码的索引零不是初始物理源节点：初始源节点是代码索引一。

### 1.2 系统核：主线选择预定义 SIAM

主线选择两个杂质模，系统 Hamiltonian 为

$$
A=\mu(n_\uparrow+n_\downarrow)+U n_\uparrow n_\downarrow,
\qquad U_h=e^{-ihA}.
\tag{1}
$$

这里代码的 $\mu$ 就是杂质能级系数，不在此另作化学势减法。令

$$
u=e^{-ih\mu},\quad X_\sigma=\bar x_\sigma y_\sigma,\quad
Y=X_\uparrow X_\downarrow
=\bar x_\uparrow\bar x_\downarrow y_\downarrow y_\uparrow.
$$

由四个 Fock 态的传播权重直接得到

$$
k_h(\bar x,y)=\langle\bar x|U_h|y\rangle
=1+u(X_\uparrow+X_\downarrow)+u^2e^{-ihU}Y
=\exp\{u(X_\uparrow+X_\downarrow)+bY\},
$$

$$
b=u^2(e^{-ihU}-1).
\tag{2}
$$

最后一个等号来自幂零性：指数中两个二次项的乘积贡献 $u^2Y$，再加 $bY$ 正好给出双占据权重。对应代码是 [siam.jl:183](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/predefined/siam.jl:183)、反分支 [siam.jl:213](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/predefined/siam.jl:213) 和系数 [siam.jl:290](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/predefined/siam.jl:290)。

通用 `ImpurityHamiltonian` 的短时核通常不是式 (2) 那样的精确传播核，另在第 10 节分析。主定理不混用这两个接口。

### 1.3 IF：完整 order=1 静态构造

[correlationfunction.jl:13](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/correlationfunction.jl:13) 调用 QuAPI 的 `Δt(bath,N=N,t=Nh)`。[fermionic.jl:92](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/realtime/fermionic.jl:92) 用步长 $h=t/N$，生成四个大小均为 $(N+1)\times(N+1)$ 的矩阵。

[partialif/realtime.jl:6](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/influencefunctional/hybridization/partialif/realtime.jl:6) 遍历 `i=1:lattice.k` 的每一行，正、反分支都包含。其 `partialif_hybrid` 收集该行所有分支的列；同文件 [realtime.jl:57](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/influencefunctional/hybridization/partialif/realtime.jl:57) 的逐项版本明确构造每个双线性指数。由于这些因子都是 Grassmann 偶元，它们相乘等于总二次型的指数。

因此需要证明的 IF 代码目标是

$$
\mathcal I_h=\exp\left[
\sum_{i,j=0}^N\sum_{\zeta,\zeta'=\pm}\sum_\alpha
\bar a_{i\alpha}^{\zeta}\Delta_{ij,\alpha}^{\zeta\zeta'}
a_{j\alpha}^{\zeta'}\right].
\tag{3}
$$

主线对应代码按 band 分别构造 IF 的接口，故假设不同杂质模的浴二点相关为对角分块。数学上允许非对角轨道相关，但不能把这个推广直接说成当前逐 band API 已经实现。

### 1.4 初态、测度、边界和观测量

默认系统初态为真空，其相干态核是常数一；见 [models.jl:23](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/models.jl:23)。更一般的粒子数守恒初态核记为

$$
r_0=\langle\bar a_0^+|\rho_S|a_0^-\rangle.
\tag{4}
$$

[sysinitstate2.jl:113](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/sysinitstate2.jl:113) 把 Fock 矩阵元放在初始正分支横线变量与初始反分支非横线变量上，正好是式 (4) 的位置。其转换函数会过滤绝对值小于默认阈值的系数，也会跳过 bra、ket 粒子数不相等的矩阵元；因此代码对应性不覆盖任意宇称偶初态。本文主定理取真空，或经此接口精确保留的粒子数守恒初态。若要对应物理归一化的分子，选 `normalize=false` 并保留初态原系数；若只改了整个初态的统一标量，归一化比值不变，但未归一化分子会相应缩放。

代码将各对 Grassmann 变量的 Gaussian 测度包含在积分操作中。本文写成

$$
d\mu(a)=\prod_\alpha d\bar a_\alpha\,da_\alpha\,
e^{-\sum_\alpha\bar a_\alpha a_\alpha}.
\tag{5}
$$

在相邻排序下，单对多项式 $c_0+c_a a+c_{\bar a}\bar a+c_2a\bar a$ 的积分为 $c_0+c_2$。这正是相邻两站点收缩时保留空、双占据两种贡献的规则；见 [integrate1.jl:39](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/integration/acintegrate/integrate1.jl:39) 和 [grassmannmps/util.jl:65](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/grassmannmps/util.jl:65)。不同张量腿方向的符号由分次收缩处理，不能另外重复加入一遍测度。

末端无观测量时，代码 [boundarycondition.jl:45](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/gvconnections/boundarycondition.jl:45) 构造

$$
T_I=\exp\left[-\sum_\alpha\bar\chi_\alpha a_{N\alpha}^+
+\sum_\alpha\bar a_{N\alpha}^-\chi_\alpha\right].
\tag{6}
$$

它的来源及负号在第 5 节从迹公式推出。完整系统因子包括这个边界、初态核以及全部系统链接；并非只指 `sysdynamics` 单独返回的对象。

对于一般偶杂质观测量，本文通过替换末端矩阵元定义其精确离散分子。package 并没有一个统一的“任意矩阵观测量”入口；对占据等已有入口，第 8 节进一步核对实际插入。本文的物理时刻 $t=Nh$ 对应代码末节点 `k=N+1`，不能直接把代码整数索引当成物理时间。

## 2. 从真实 Hamiltonian 和真实物理目标出发

现在固定有限个浴模，所有算符有界。完整 CAR 包括杂质奇算符与浴奇算符之间的反对易关系。令

$$
H=A+B,\qquad
B=H_B+\sum_\alpha(d_\alpha^\dagger F_\alpha+F_\alpha^\dagger d_\alpha),
$$

$$
H_B=\sum_{\alpha k}\epsilon_{\alpha k}b_{\alpha k}^\dagger b_{\alpha k},
\qquad F_\alpha=\sum_k v_{\alpha k}b_{\alpha k}.
\tag{7}
$$

$A$ 主线取式 (1)。不同 $\alpha$ 使用独立浴模，初态为

$$
\rho_0=\rho_S\otimes\rho_B,\qquad
\operatorname{Tr}\rho_S=\operatorname{Tr}\rho_B=1.
\tag{8}
$$

两部分初态均为宇称偶。浴态平稳、正规 Gaussian，满足 Wick 定理且无反常配对；为直接对应 QuAPI，取自由浴的热平衡态，各模占据为 $f_{\alpha k}$。杂质初态按第 1.4 节选择。观测量 $O$ 作用于杂质且宇称偶，其物理目标是

$$
W(t)=e^{-itH},\qquad
N_O(t)=\operatorname{Tr}_{SB}[Oe^{-itH}\rho_0e^{+itH}],
\qquad N_I(t)=1.
\tag{9}
$$

真实物理期望值就是 $N_O(t)$。离散算法的分母一般不严格等于一，所以还需要证明它的归一化比值收敛。

我们的策略是先找出一个来自式 (9)、具有算符范数误差界的近似传播子，再证明第 1 节的多项式恰好表示它。为了构造这个传播子，引入下一节的辅助映射。

## 3. 单独建立重建映射 $Q$

### 3.1 为什么引入它，以及固定什么约定

代码 IF 的同一节点同时带 $\bar a_j,a_j$，而系统核连接不同节点。普通的两传播子相干态核拼接不足以推出这种结构。为重建当前算法，我们选择一个对角相干态框架映射

$$
Q[f]=\int d\mu(a)|a\rangle f(\bar a,a)\langle\bar a|.
\tag{10}
$$

这是辅助选择：package 没有实际计算这个算符积分，原始 Hamiltonian 也不唯一规定必须作这个选择。后面将证明：这个选择一阶逼近真实耦合传播，而且精确生成所核对的代码目标。

固定单模 Fock 态 $|0\rangle,|1\rangle=d^\dagger|0\rangle$，使用未归一化相干态

$$
|a\rangle=e^{-a d^\dagger}|0\rangle,
\qquad \langle\bar a|=\langle0|e^{-d\bar a},
\qquad \langle\bar x|y\rangle=e^{\bar x y}.
\tag{11}
$$

所有 Grassmann 变量独立，横线不是数值复共轭。它们与所有费米奇算符及奇基矢反对易。张量积使用 Koszul 规则：交换两个齐次对象的符号是其宇称乘积的负一次幂。相干态 ket、bra 本身总宇称为偶，例如 $-a|1\rangle$ 是奇系数乘奇基矢。

单模规范顺序取 $a,\bar a$，规定 $\int d\bar a\,da\,a\bar a=1$；多模规范单项式先按模排列为 $a_1,\bar a_1,a_2,\bar a_2,\ldots$。物理 Fock 基固定为 $(d_1^\dagger)^{n_1}\cdots(d_m^\dagger)^{n_m}|0\rangle$，所有杂质模排在浴模之前，浴模内部也固定一个顺序。

完整轮廓多项式的规范生成元顺序取：先是按模排列的 $\chi_\alpha,\bar\chi_\alpha$；随后节点 $j$ 从零递增，每个节点先正分支、后反分支，每个分支内再按模排列 $a_{j\alpha}^\zeta,\bar a_{j\alpha}^\zeta$。浴系数统一排在全部 Grassmann 变量右侧。每个积分变量对为偶，变量对之间的积分可交换。代码默认 MPS 次序可由 [realtime.jl:153](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/lattices/realtime.jl:153) 的 `index` 函数读出，它与这里的规范排列不必相同；转换单项式时必须同步使用 Grassmann 置换符号，不能只改节点名字。

### 3.2 单模逐项计算，不借用传播子拼接结论

取独立于 $a,\bar a$ 的外部变量 $\bar x,y$。按式 (10)—(11)，

$$
\langle\bar x|Q[f]|y\rangle
=\int d\bar a\,da\,(1+a\bar a)
(1+\bar x a)f(\bar a,a)(1+\bar a y).
\tag{12}
$$

逐项把积分变量移到顺序 $a\bar a$，抽取其系数：

| $f$ | 式 (12) 中 $a\bar a$ 的系数 | 算符核的识别 |
|---|---|---|
| $1$ | $1+\bar x y$ | $\langle\bar x|I|y\rangle$ |
| $a$ | $y$ | $\langle\bar x|d|y\rangle$ |
| $\bar a$ | $\bar x$ | $\langle\bar x|d^\dagger|y\rangle$ |
| $a\bar a$ | $1$ | $\langle\bar x|(I-n)|y\rangle$ |

例如第二行只有 $a\bar a y$ 贡献积分；第三行的贡献是 $\bar x a\bar a$；第四行已有全部积分变量，其余含它们的项都因平方为零而消失。第一行由测度的 $a\bar a$ 项和两个重叠因子的乘积分别给出 $1$ 与 $\bar x y$。因此

$$
Q[1]=I,\quad Q[a]=d,\quad Q[\bar a]=d^\dagger,\quad
Q[a\bar a]=I-n,\quad Q[\bar a a]=n-I.
\tag{13}
$$

为什么核相同就能确定算符？相干态在 Fock 基中的展开把每个矩阵元乘上一项独立的外部 Grassmann 单项式。零核的每个系数都为零，故所有 Fock 矩阵元为零。多模时 bra 占据子集和 ket 占据子集仍分别标记独立的单项式，论证不变。因此式 (13) 是算符等式，不只是一组特定态平均值。

式 (13) 同时验证 $Q[1]=I$ 的恒等分解。对普通宇称偶算符逐项计算还给出

$$
\operatorname{Tr}_S X=\int d\mu(\chi)\langle-\bar\chi|X|\chi\rangle.
\tag{14}
$$

单模对角矩阵元的核为 $X_{00}-\bar\chi\chi X_{11}$，乘测度后积分得到 $X_{00}+X_{11}$。若去掉负号则得到宇称加权迹，故这个负号不能省略。

### 3.3 含浴奇系数时的符号

把一个总宇称偶的单模多项式唯一写为

$$
f=c_0+a c_a+\bar a c_{\bar a}+a\bar a c_2,
\tag{15}
$$

浴系数统一放在右侧。$c_0,c_2$ 为偶，$c_a,c_{\bar a}$ 为奇。将式 (15) 原样代入式 (12)，保留系数与外部奇变量的反对易关系，结果是

$$
\langle\bar x|Q[f]|y\rangle
=(1+\bar x y)c_0+y c_a+\bar x c_{\bar a}+c_2,
$$

$$
Q[f]=c_0+d c_a+d^\dagger c_{\bar a}+(I-n)c_2.
\tag{16}
$$

比如相关项为 $a c_a\bar a y$：先移动 $c_a$ 经过 $\bar a$，再经过 $y$，两个负号抵消，成为 $a\bar a y c_a$。结果是 $y c_a$，不能把中间的浴奇系数当成普通数。

### 3.4 多模推广：给出完整单项式规则

对规范单项式 $M=\prod_{\alpha=1}^m m_\alpha$，每个 $m_\alpha$ 取 $1,a_\alpha,\bar a_\alpha,a_\alpha\bar a_\alpha$ 之一。令其单模像依次为 $I,d_\alpha,d_\alpha^\dagger,I-n_\alpha$，记作 $q_\alpha(m_\alpha)$。则

$$
Q[M c]=q_1(m_1)q_2(m_2)\cdots q_m(m_m)c,
\tag{17}
$$

仍将浴系数 $c$ 放在最右。

下面说明这个规则为什么成立。按模作分次张量积，相干态及投影框架都是单模对象的张量积；每对测度均为偶。先积分第一个模，式 (12) 将其四种可能因子替换成式 (13) 的算符，剩余模及浴系数作为固定的分次系数保留。随后逐模重复。同一次重排中，若一个奇系数必须越过某个奇模因子，积分前产生的负号在替换后完全相同，因为式 (13) 保持各基因子的宇称。这样归纳得到式 (17)。它覆盖全部单项式，再由线性覆盖全部多项式。

若希望用普通矩阵而非分次张量积复核，令 $P_S=(-1)^{\sum n_\alpha}$。对一个总偶项 $Mc$，浴系数宇称等于 $M$ 的次数模二，记为 $p$。式 (17) 的普通张量积矩阵正是

$$
\bigl[q_1(m_1)\cdots q_m(m_m)P_S^p\bigr]\otimes c.
\tag{18}
$$

这里杂质模算符自身采用固定 Fock 排序的 Jordan–Wigner 矩阵。式 (18) 的宇称串保证浴奇算符与全部杂质奇算符反对易；省略它会给出错误的多模重建。

后面需要的特例为

$$
\begin{aligned}
Q[a_\alpha\bar a_\beta]&=d_\alpha d_\beta^\dagger,\qquad
Q[\bar a_\alpha a_\beta]=d_\alpha^\dagger d_\beta-\delta_{\alpha\beta}I,\\
Q[a_\alpha a_\beta]&=d_\alpha d_\beta,\qquad
Q[\bar a_\alpha\bar a_\beta]=d_\alpha^\dagger d_\beta^\dagger.
\end{aligned}
\tag{19}
$$

伴随约定为 $a_\alpha^\dagger=\bar a_\alpha$，并反转整个乘积次序。单模表逐项满足伴随相容；多模伴随后将模因子和浴系数移回规范顺序时，两边因宇称相同而得到相同符号。因此

$$
Q[f^\dagger]=Q[f]^\dagger.
\tag{20}
$$

例如 $(a c)^\dagger=c^\dagger\bar a=-\bar a c^\dagger$ 对应 $(dc)^\dagger=c^\dagger d^\dagger=-d^\dagger c^\dagger$。这一步将在构造反分支传播子时使用。

$Q$ 是线性映射，但不保乘法：$Q[\bar a a]=n-I$，而 $Q[\bar a]Q[a]=n$。这一差别正是源传播重建的局域误差来源。

## 4. 先完成真实传播到 order=1 近似传播的误差控制

### 4.1 精确定义源单元，再作 Taylor 分析

为使 IF 最终在每个节点同时含两个源，定义

$$
L(\bar a,a)=H_B+\sum_\alpha(\bar a_\alpha F_\alpha+F_\alpha^\dagger a_\alpha),
\qquad E_w=e^{-iwL},\qquad R_w=Q[E_w].
\tag{21}
$$

式 (21) 都是精确定义，尚未把真实传播子替换成它们。$L$ 是带外部 Grassmann 源的浴算符，已经不含杂质算符；它不是原来的 $B$。由 $F_\alpha^\dagger a_\alpha=-a_\alpha F_\alpha^\dagger$ 及式 (16)，

$$
Q[L]=H_B+\sum_\alpha(d_\alpha^\dagger F_\alpha+F_\alpha^\dagger d_\alpha)=B.
\tag{22}
$$

固定有限浴，外代数也有限维。可用“规范单项式各浴系数范数之和”作为次乘性范数；式 (17)—(18) 说明 $Q$ 对此范数有界。故指数的 Taylor 余项可控制，得到

$$
R_w=I-iwB-\frac{w^2}{2}Q[L^2]+O(w^3).
\tag{23}
$$

这没有把 $R_w$ 等同于 $I-iwB$。

### 4.2 重新算出二次差额

含 $H_B$ 的各项在重建后与 $B^2$ 的相应项相同，因为 $H_B$ 为偶浴算符。纯耦合项逐类比较：

| 源中相乘的两个顶点 | 提出源后的表达式 | 重建与真实算符乘积的差额 |
|---|---|---|
| $(\bar a_\alpha F_\alpha)(\bar a_\beta F_\beta)$ | $-\bar a_\alpha\bar a_\beta F_\alpha F_\beta$ | $0$ |
| $(F_\alpha^\dagger a_\alpha)(F_\beta^\dagger a_\beta)$ | $-a_\alpha a_\beta F_\alpha^\dagger F_\beta^\dagger$ | $0$ |
| $(F_\alpha^\dagger a_\alpha)(\bar a_\beta F_\beta)$ | $a_\alpha\bar a_\beta F_\alpha^\dagger F_\beta$ | $0$ |
| $(\bar a_\alpha F_\alpha)(F_\beta^\dagger a_\beta)$ | $\bar a_\alpha a_\beta F_\alpha F_\beta^\dagger$ | $-\delta_{\alpha\beta}F_\alpha F_\beta^\dagger$ |

最后一行使用式 (19)，其他三行也分别由式 (19) 验证。故

$$
Q[L^2]=B^2-G,\qquad
G=\sum_\alpha F_\alpha F_\alpha^\dagger\geq0,
$$

$$
\boxed{R_w=e^{-iwB}+\frac{w^2}{2}G+O(w^3).}
\tag{24}
$$

现在才实际选择用 $R_h$ 替换真实的 $V_h=e^{-ihB}$。对充分小的 $h$，式 (24) 给出常数 $c$，使

$$
\|R_h-V_h\|\leq ch^2,\qquad \|R_h\|\leq1+ch^2.
\tag{25}
$$

常数可依赖固定有限浴与模型参数，但不依赖 $N=t/h$。$R_h$ 一般不幺正。

### 4.3 真实演化 → Lie 分裂 → 额外端点 → 逐个源重建

这部分先证明近似是可控的，再进入其路径表示。定义真实的系统、耦合块短时传播子 $U_h=e^{-ihA}$、$V_h=e^{-ihB}$。

第一步，使用固定顺序的 Lie 乘积

$$
S_h=(V_hU_h)^N.
$$

两次 Duhamel 积分把单步差额写成带 $[A,B]$ 的双积分；其余因子都是幺正的。因此

$$
\|e^{-ih(A+B)}-V_hU_h\|\leq\frac{h^2}{2}\|[A,B]\|,
\quad
\|W(t)-S_h\|\leq\frac{th}{2}\|[A,B]\|.
\tag{26}
$$

第二步，为匹配代码的 $N+1$ 个整宽节点，在最右侧加入一个真实的 $B$ 传播子：

$$
T_h=S_hV_h=(V_hU_h)^NV_h.
$$

它从右到左先作用于初态，因此这里明确选定的是额外的初始端点块。其代价是

$$
\|T_h-S_h\|=\|V_h-I\|\leq h\|B\|.
\tag{27}
$$

不能跳过式 (27)。$T_h$ 并非仍只含 $N$ 个 $B$ 传播块。

第三步，把 $T_h$ 中全部 $N+1$ 个 $V_h$ 逐个换成 $R_h$，定义

$$
\boxed{P_h=(R_hU_h)^NR_h
=R_hU_hR_hU_h\cdots U_hR_h.}
\tag{28}
$$

乘积差的伸缩求和中，每一项含一个 $R_h-V_h$，其余因子范数由式 (25) 控制，所以

$$
\|P_h-T_h\|\leq(N+1)ch^2(1+ch^2)^N
\leq ch(t+h)e^{cth}.
\tag{29}
$$

合并三步，得到完整的误差链

$$
\boxed{
W(t)\ \xrightarrow{\ O(h)\ }\ S_h
\ \xrightarrow{\ O(h)\ }\ T_h
\ \xrightarrow{\ O(h)\ }\ P_h,
\qquad \|P_h-W(t)\|\leq C_t h.}
\tag{30}
$$

各箭头分别有式 (26)、(27)、(29) 的独立上界。系统传播间隔总长仍是 $Nh=t$，源块总长是 $(N+1)h=t+h$。二者并不矛盾：额外的一个 $B$ 块本身只造成 $O(h)$ 差额。也不能把这个混合乘积直接解释为真实 Hamiltonian 演化到 $t+h$，因为系统块没有增加。

### 4.4 将误差传到密度矩阵、分子和归一化结果

此时已证明 $P_h$ 逼近真实传播子，可以定义其输出而不循环论证。令

$$
\rho_h=P_h\rho_0P_h^\dagger,\quad
\rho(t)=W(t)\rho_0W(t)^\dagger,\quad
N_{O,h}=\operatorname{Tr}_{SB}(O\rho_h),\quad Z_h=\operatorname{Tr}_{SB}\rho_h.
\tag{31}
$$

用 $\|\rho_0\|_1=1$ 和 $\|W(t)\|=1$，

$$
\|\rho_h-\rho(t)\|_1
\leq\|P_h-W(t)\|(\|P_h\|+1)=O(h).
\tag{32}
$$

约化密度矩阵也满足这个界，因为部分迹对自伴差额的迹范数不增。于是

$$
|N_{O,h}-N_O(t)|\leq\|O\|\,\|\rho_h-\rho(t)\|_1=O(h),
\qquad |Z_h-1|=O(h).
\tag{33}
$$

充分小步长时 $Z_h\geq1/2$，且

$$
\left|\frac{N_{O,h}}{Z_h}-N_O(t)\right|
\leq\frac{|N_{O,h}-N_O(t)|+|N_O(t)|\,|Z_h-1|}{Z_h}
=O(h).
\tag{34}
$$

这里仅证明一阶上界，尚未证明一般不具有二阶；后者在第 9 节用当前整宽端点版本单独分析。

## 5. 将已受控的近似算符精确展开为共同 Grassmann 积分

### 5.1 一个系统间隔：完整保留四套节点变量

先取 $N=1$，为每个节点和分支分配独立变量。记

$$
E_j^+=e^{-ihL(\bar a_j^+,a_j^+)},\qquad
E_j^-=e^{+ihL(\bar a_j^-,a_j^-)},\qquad j=0,1.
$$

由 $Q$ 的定义及伴随相容性，直接代入而不作任何节点识别，得到

$$
P_h=\int d\mu(a_1^+)d\mu(a_0^+)\,
|a_1^+\rangle E_1^+
\underbrace{\langle\bar a_1^+|U_h|a_0^+\rangle}_{k_0^+}
E_0^+\langle\bar a_0^+|,
$$

$$
P_h^\dagger=\int d\mu(a_0^-)d\mu(a_1^-)\,
|a_0^-\rangle E_0^-
\underbrace{\langle\bar a_0^-|U_h^\dagger|a_1^-\rangle}_{k_0^-}
E_1^-\langle\bar a_1^-|.
\tag{35}
$$

把这两式和因子化初态放进同一个乘积，

$$
P_h\rho_0P_h^\dagger
=\int d\nu_1\,
|a_1^+\rangle E_1^+k_0^+E_0^+
\langle\bar a_0^+|\rho_S|a_0^-\rangle\rho_B
E_0^-k_0^-E_1^-\langle\bar a_1^-|,
$$

$$
d\nu_1=d\mu(a_0^+)d\mu(a_1^+)d\mu(a_0^-)d\mu(a_1^-).
\tag{36}
$$

系统传播、初态均为偶算符，故 $k_0^+,k_0^-,r_0$ 都是纯 Grassmann 偶多项式，可以移到浴算符乘积之外。$E_j$ 虽然总宇称偶，却共享同一套浴算符，因此不同 $E_j$ 一般不对易；也不能把它们任意穿过 $\rho_B$。

做杂质迹时，式 (14) 引入一套独立的辅助变量。对偶观测量 $O$，定义

$$
T_O=\langle-\bar\chi|O|a_1^+\rangle
\langle\bar a_1^-|\chi\rangle.
\tag{37}
$$

它只负责完成 Fock 迹，未增加物理传播时间。代入 $O=I$，两个重叠给出式 (6) 的负号和正号。于是式 (36) 严格变成

$$
N_{O,h}=\int d\mu(\chi)d\nu_1\,
\underbrace{T_O r_0 k_0^+k_0^-}_{K_{O,h}}
\underbrace{\operatorname{Tr}_B[E_1^+E_0^+\rho_B E_0^-E_1^-]}_{\mathcal I_h}.
\tag{38}
$$

式 (35)—(38) 都是受控近似传播子 $P_h$ 的恒等展开。

### 5.2 一般 $N$、完整定义和共享变量的来源

对 $j=0,\ldots,N-1$，定义

$$
k_j^+=\langle\bar a_{j+1}^+|U_h|a_j^+\rangle,\qquad
k_j^-=\langle\bar a_j^-|U_h^\dagger|a_{j+1}^-\rangle.
$$

重复式 (35) 的替换，保持所有节点变量独立，得到

$$
\begin{aligned}
K_{O,h}&=T_O\,r_0\prod_{j=0}^{N-1}k_j^+k_j^-,\\
T_O&=\langle-\bar\chi|O|a_N^+\rangle\langle\bar a_N^-|\chi\rangle,\\
\mathcal I_h&=\operatorname{Tr}_B[E_N^+\cdots E_0^+\rho_B E_0^-\cdots E_N^-],\\
d\nu_h&=d\mu(\chi)\prod_{j=0}^N d\mu(a_j^+)d\mu(a_j^-).
\end{aligned}
\tag{39}
$$

这些纯 Grassmann 偶因子之间可交换，但所有浴单元在式 (39) 中的次序固定。最终精确得到

$$
\boxed{N_{O,h}=\int d\nu_h\,K_{O,h}\mathcal I_h,\qquad
Z_h=\int d\nu_h\,K_{I,h}\mathcal I_h.}
\tag{40}
$$

看一个内部正分支节点 $j$：$E_j^+$ 含 $a_j^+,\bar a_j^+$；$k_j^+$ 使用 $a_j^+$，$k_{j-1}^+$ 使用 $\bar a_j^+$。所以 $K$ 与 IF 同时包含这一对变量。反分支同理；端点则由初态和迹边界补齐。

这种交集来自 $|a_j\rangle E_j(\bar a_j,a_j)\langle\bar a_j|$ 的重建结构。它没有把独立节点强行等同，也不是直接套用两个精确传播子相干态核的互补接缝公式。将它作为真实传播子的替代所需代价，已经由式 (24)—(34) 单独证明。

## 6. 从 E 的乘积浴迹直接计算 IF

### 6.1 连续时间从有限指数的精确 Dyson 展开出现

直接从式 (39) 的浴迹开始。因为每个源单元的宽度都是 $h$，为浴演化引入连续区间

$$
C_j=[jh,(j+1)h],\quad j=0,\ldots,N,\qquad T=(N+1)h=t+h.
\tag{41}
$$

这些是按 $E_N\cdots E_0$ 的先后排列定义的浴源时钟。节点的系统标记时间仍是 $jh$；区间的左端被标为这个节点，不意味着源只有点值作用。平稳浴下把所有区间共同平移不改变最终系数。

在 $s\in C_j$ 上定义

$$
V_j^\zeta(s)=\sum_\alpha
[\bar a_{j\alpha}^\zeta F_\alpha(s)+F_\alpha^\dagger(s)a_{j\alpha}^\zeta],
\quad F_\alpha(s)=e^{isH_B}F_\alpha e^{-isH_B}.
\tag{42}
$$

将每个有限时间指数 $E_j$ 精确转到相互作用绘景，再合并相邻区间，得到

$$
E_N^+\cdots E_0^+=e^{-iTH_B}W_+,
\quad W_+=\mathcal T\exp\left[-i\sum_j\int_{C_j}V_j^+(s)\,ds\right],
$$

$$
E_0^-\cdots E_N^-=\widetilde W_-e^{+iTH_B},
\quad\widetilde W_-=\widetilde{\mathcal T}\exp\left[+i\sum_j\int_{C_j}V_j^-(s)\,ds\right].
\tag{43}
$$

$\mathcal T$ 把较晚的正分支浴时间放左侧，$\widetilde{\mathcal T}$ 把较早的负分支浴时间放左侧。这里的时间排序作用于总偶顶点，无需另加顶点交换的费米负号；内部奇场的 Wick 符号仍须保留。

式 (43) 来自有限指数的 Dyson 恒等式，不是把真实杂质路径额外近似成连续函数。源在单元内固定这一选择已经包含于式 (21)。外侧自由浴演化是偶算符；逐个 Grassmann 系数取普通浴矩阵迹，可消去外侧相似变换，因而

$$
\mathcal I_h=\operatorname{Tr}_B[W_+\rho_B\widetilde W_-].
\tag{44}
$$

这只使用了偶自由浴算符的迹不变性，没有把带奇源系数的整个表达式循环搬动。

### 6.2 先算二点项，固定四个分支核和符号

省略对角轨道标记，定义无额外 $i$ 因子的两个关联函数

$$
C^>(s,s')=\operatorname{Tr}_B[\rho_B F(s)F^\dagger(s')],\qquad
C^<(s,s')=\operatorname{Tr}_B[\rho_B F^\dagger(s')F(s)].
\tag{45}
$$

有限浴中它们分别为

$$
C^>(s,s')=\sum_k|v_k|^2(1-f_k)e^{-i\epsilon_k(s-s')},\quad
C^<(s,s')=\sum_k|v_k|^2f_ke^{-i\epsilon_k(s-s')}.
\tag{46}
$$

正分支的二阶 Dyson 项在 $s>s'$ 上带 $(-i)^2=-1$。两个非零混合顶点分别为

$$
(\bar a(s)F(s))(F^\dagger(s')a(s'))
=\bar a(s)a(s')F(s)F^\dagger(s'),
$$

$$
(F^\dagger(s)a(s))(\bar a(s')F(s'))
=-\bar a(s')a(s)F^\dagger(s)F(s').
$$

所以正分支下三角的系数是 $-C^>$，上三角的系数是 $+C^<$。

负分支在同一区域 $s>s'$ 使用反时间次序 $V^-(s')V^-(s)$，仍有 $(+i)^2=-1$。其两个混合项给出

$$
-\bar a^-(s')a^-(s)C^>(s',s)
+\bar a^-(s)a^-(s')C^<(s,s').
$$

因此负分支下三角是 $+C^<$，上三角是 $-C^>$。

跨分支项没有时间三角限制，且 $(-i)(+i)=1$。保留密度矩阵的位置，展开

$$
\operatorname{Tr}_B[
(\bar a^+(s)F(s)+F^\dagger(s)a^+(s))\rho_B
(\bar a^-(s')F(s')+F^\dagger(s')a^-(s'))].
$$

先把外部源按分次规则移到前面，再对剩下的纯浴矩阵取迹，得到

$$
\bar a^+(s)a^-(s')\operatorname{Tr}_B[F(s)\rho_BF^\dagger(s')]
=\bar a^+(s)a^-(s')C^<(s,s'),
$$

$$
-\bar a^-(s')a^+(s)\operatorname{Tr}_B[F^\dagger(s)\rho_BF(s')]
=-\bar a^-(s')a^+(s)C^>(s',s).
\tag{47}
$$

最后两个等号只循环了已去除源的普通浴算符。跨分支第二项的负号来自源排序，不能漏掉。

因此全部二点源项可写为

$$
D=\sum_{\zeta\zeta'}\sum_{ij}\int_{C_i}ds\int_{C_j}ds'\,
\bar a_i^\zeta\Xi^{\zeta\zeta'}(s,s')a_j^{\zeta'},
$$

$$
\begin{aligned}
\Xi^{++}(s,s')&=-\theta(s-s')C^>(s,s')+\theta(s'-s)C^<(s,s'),\\
\Xi^{--}(s,s')&=\theta(s-s')C^<(s,s')-\theta(s'-s)C^>(s,s'),\\
\Xi^{+-}(s,s')&=C^<(s,s'),\\
\Xi^{-+}(s,s')&=-C^>(s,s').
\end{aligned}
\tag{48}
$$

对角轨道只需分别相加。有限浴的核有界，时间相等的直线为零测集，阶跃函数在零点的取值不影响这些积分；同一单元的上下三角却都有非零面积，必须都计算。

### 6.3 为什么全部高阶项恰好等于 $\exp(D)$

这里必须证明全阶重求和，不能在二点项算完后直接宣称成立。

在式 (44) 的全部耦合顶点前乘一个标量 $\lambda$，按 Dyson 级数展开。每个顶点线性含一个浴奇场和一个外部奇源。固定正、负分支的顶点数后，源提出到前方，剩余的是形如 $\operatorname{Tr}_B[X_1\cdots X_p\rho_B Y_1\cdots Y_q]$ 的纯浴奇场迹。此时可以先把普通浴乘积循环成 $\operatorname{Tr}_B[\rho_B Y_1\cdots Y_qX_1\cdots X_p]$，再使用 Gaussian 浴的 Wick 定理。这个操作在每个系数上进行，不是循环带 Grassmann 系数的整个浴迹。

奇数顶点的浴迹为零。对 $2r$ 个顶点，Wick 定理将浴场期望写为所有完全二点配对之和，每项含相应的费米置换符号。

先核实密度矩阵夹心形式本身不破坏这个配对规则。对固定的 $p+q=2r$ 个浴奇场，取原标签次序为全部 $X$ 后接全部 $Y$。把浴迹循环成以 $\rho_B$ 开头后，场次序变为全部 $Y$ 后接全部 $X$。这次块置换的奇偶性是 $(-1)^{pq}$。将普通 Wick 配对重新用原标签排序时，每个跨块配对的协方差相对原标签的夹心二点迹多一个负号：原标签夹心二点迹是 $\operatorname{Tr}(X\rho_BY)=\langle YX\rangle_B$，而循环后协方差在恢复原标签时因反对称性成为 $-\langle YX\rangle_B$。若一个完全配对有 $c$ 个跨块配对，则 $p-c$、$q-c$ 都为偶数，故 $c\equiv p\equiv pq\pmod2$。所以块置换符号与全部跨块符号相乘为 $(-1)^{pq+c}=1$。由此，原夹心迹确实按它自身的夹心二点迹满足 Wick 配对展开；没有暗中使用错误的 Grassmann 迹循环性。

对每一种配对，把对应源也排成相同的成对次序。浴场排序符号与源排序符号一起，恰好把该项变成相应两点“源乘收缩”因子的乘积。不同的两点源因子为偶元，故彼此可交换。四种分支二点源因子的具体值就是式 (48)。

可以用四点的基本结构检查这一步：浴 Wick 符号是 $C_{12}C_{34}-C_{13}C_{24}+C_{14}C_{23}$；把源分组成 $(12)(34)$、$(13)(24)$、$(14)(23)$ 时，源的置换正好生成同样的三种代数次序。于是应当相加的是三个两点源因子的乘积；不能把浴配对负号与源置换负号只保留其中一组。

接下来处理组合系数。正、负分支 Dyson 因子分别带 $1/p!$、$1/q!$；把分支、时间、轨道和两种源类型一起视为一个顶点标签求和后，它们等价于总共 $2r$ 个标签的 $1/(2r)!$，乘以带分支标签的排序规则。两点排序后的源因子关于两个标签的互换是对称的，因为被排序的是总偶顶点。

完全配对的数目是 $(2r)!/(2^r r!)$。在每一配对中，标签为独立的积分哑变量，因而对所有标签积分后的结果相同：每一对给出两点标签积分 $2D$。所以 $2r$ 阶项严格为

$$
\frac{1}{(2r)!}\frac{(2r)!}{2^r r!}(2D)^r
=\frac{D^r}{r!}.
\tag{49}
$$

这也是 Wick 配对组合数在 Grassmann 源代数中的通常指数公式；上述推导说明它适用于式 (44) 的密度矩阵夹心次序，并保留了跨分支符号。零顶点项等于 $\operatorname{Tr}\rho_B=1$，所以没有额外的常数行列式。

从而逐阶得到

$$
\mathcal I_h(\lambda)=\sum_{r\geq0}\frac{\lambda^{2r}D^r}{r!}
=e^{\lambda^2D},\qquad \boxed{\mathcal I_h=e^D.}
\tag{50}
$$

这里独立源只有 $4m(N+1)$ 个，超过这个次数的项自动为零。因此式 (49)—(50) 实际是有限 Grassmann 多项式的逐系数恒等式，没有把杂化强度截断在二阶，也没有依赖难以控制的无穷阶极限交换。

## 7. 从这些精确时间积分逐项得到 QuAPI 系数

### 7.1 所有单元都是整宽 $h$

将式 (48) 的二点项按式 (41) 分组，定义

$$
\boxed{\Delta_{ij,\alpha}^{\zeta\zeta'}
=\int_{jh}^{(j+1)h}\!ds'\int_{ih}^{(i+1)h}\!ds\,
\Xi_\alpha^{\zeta\zeta'}(s,s').}
\tag{51}
$$

于是由式 (50) 精确得到

$$
\boxed{\mathcal I_h=
\exp\left[\sum_{i,j=0}^N\sum_{\zeta\zeta',\alpha}
\bar a_{i\alpha}^{\zeta}\Delta_{ij,\alpha}^{\zeta\zeta'}
a_{j\alpha}^{\zeta'}\right].}
\tag{52}
$$

式 (51) 对连续核作矩形或三角区域积分，就是当前 QuAPI 的离散系数。它来自有限指数的精确 Dyson 展开和重求和；没有先构造一个不同的离散浴 Gaussian 矩阵再取 Schur 补来冒充此系数。

### 7.2 单浴能级的矩形和三角积分

对单个能量 $\epsilon$，从式 (46) 提出耦合平方及占据权重。不同单元的时间差符号固定，所需矩形积分是

$$
W_{ij}(\epsilon)=\int_{C_i}ds\int_{C_j}ds'\,e^{-i\epsilon(s-s')}
=e^{-i\epsilon(i-j)h}\frac{2(1-\cos\epsilon h)}{\epsilon^2}.
\tag{53}
$$

同一单元的正时间差部分是

$$
T_h(\epsilon)=\int_0^hds\int_0^sds'\,e^{-i\epsilon(s-s')}
=\frac{1-i\epsilon h-e^{-i\epsilon h}}{\epsilon^2}.
\tag{54}
$$

负时间差部分为 $T_h^*$。当 $\epsilon\to0$，两者极限分别为 $h^2$ 与 $h^2/2$。这对应 QuAPI 的 [_gⱼₖ_r:44](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/realtime/util.jl:44) 和 [_gⱼⱼ_r:60](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/realtime/util.jl:60)；非平凡色散通过同文件的包装函数代入。

单能级四块矩阵的完整结果如下，整张表乘以 $|v|^2$。

| 分支块 | $i>j$ | $i<j$ | $i=j$ |
|---|---|---|---|
| $++$ | $-(1-f)W_{ij}$ | $fW_{ij}$ | $-(1-f)T_h+fT_h^*$ |
| $--$ | $fW_{ij}$ | $-(1-f)W_{ij}$ | $fT_h-(1-f)T_h^*$ |
| $+-$ | $fW_{ij}$ | $fW_{ij}$ | $f(T_h+T_h^*)$ |
| $-+$ | $-(1-f)W_{ij}$ | $-(1-f)W_{ij}$ | $-(1-f)(T_h+T_h^*)$ |

跨分支项没有三角排序限制，即使 $i=j$ 也必须用完整矩形。所有 $i=0$ 和 $i=N$ 的端点单元仍为整宽 $h$，对角项也用这张表，不能把端点三角面积再减半。

有限浴的 $\Delta$ 是各能级贡献之和。这已经完成有限模型中的推导。QuAPI 的谱版本把同一张表对 $J(\epsilon)$ 积分，其 $g_1=1-f$、$g_2=f$ 的定义见 [QuAPI/util.jl:4](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/util.jl:4)，四块谱求积见 [fermionic.jl:100](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/realtime/fermionic.jl:100)。这说明系数公式逐项相同；并不把 delta 离散谱误称为 `AbstractBoundedFunction` 接口可以直接输入的普通有界函数。有限浴可直接使用有限模求和形成相同的相关矩阵。

### 7.3 对角元两个存储项相加，不是丢掉另一个三角

QuAPI 的 `CorrelationMatrix` 存下三角与上三角两个数组。其 [corrmat.jl:18](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/corrmat.jl:18) 定义

$$
\Delta_{ij}=\begin{cases}
\eta_{\rm lower}[i-j+1],&i>j,\\
\eta_{\rm upper}[j-i+1],&i<j,\\
\eta_{\rm lower}[1]+\eta_{\rm upper}[1],&i=j.
\end{cases}
\tag{55}
$$

最后一行尤其重要。比如 $+-$ 的两个首项分别存 $fT_h$ 和 $fT_h^*$，相加才是完整矩形；$-+$ 同理。$++$、$--$ 则分别相加同单元中按不同占据权重计算的两个三角。代码没有把跨分支对角矩形少算一半。

由第 1.3 节的逐行构造，式 (52) 的全部系数恰好被放在同一节点、分支、band 的变量上。因此到这里才完成“浴迹等于当前 IF”的对应证明。

## 8. 将实际代码、物理观测量和主定理接起来

### 8.1 完整系统因子及无观测量积分

对主线 SIAM，式 (2) 正是 $k_j^+$，其伴随是 $k_j^-$。式 (4) 对应初态，式 (6) 对应迹边界。因此代码无观测量时构造的完整系统因子为

$$
K_h=T_I r_0\prod_{j=0}^{N-1}k_h(\bar a_{j+1}^+,a_j^+)
\langle\bar a_j^-|U_h^\dagger|a_{j+1}^-\rangle.
\tag{56}
$$

测度由积分例程执行，不再包含进式 (56)。如果把初态和边界存成独立 GMPS，也必须在积分中一起乘入。代码的 [integration.jl:37](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/integration/integration.jl:37) 对输入 GMPS 的乘积积分；使用 `ExactIntegrate` 时不在这一步引入新的积分截断。

由式 (40)、(52)、(56)，

$$
Z_h=\operatorname{Tr}[P_h\rho_0P_h^\dagger]
=\int d\nu_h\,K_h\mathcal I_h
=Z_{\rm code,h}.
\tag{57}
$$

末尾等号指已核对系数的代码离散目标，不是声称执行 Julia 后的有限精度结果绝无误差。

### 8.2 末时刻占据与实际 Gt 插入

对一般偶 $O$，式 (39) 的 $T_O$ 已经给出物理分子 $N_{O,h}$ 的精确表示。对已有占据数入口，可进一步消去这个特殊边界核。

相干态满足 $d_\alpha|a_N^+\rangle=a_{N\alpha}^+|a_N^+\rangle$、$\langle\bar a_N^-|d_\alpha^\dagger=\langle\bar a_N^-|\bar a_{N\alpha}^-$。故

$$
\operatorname{Tr}[n_\alpha P_h\rho_0P_h^\dagger]
=\operatorname{Tr}[d_\alpha P_h\rho_0P_h^\dagger d_\alpha^\dagger]
=\int d\nu_h\,a_{N\alpha}^+\bar a_{N\alpha}^-\,K_h\mathcal I_h.
\tag{58}
$$

第一个等号是完整物理算符的普通迹循环性，尚无外部 Grassmann 系数；第二个等号通过式 (35) 两端的相干态本征关系得到，保持源的次序 $a_N^+\bar a_N^-$。

[occupation.jl:8](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/occupation.jl:8) 调用正分支非横线、反分支横线的 `Gt(i,i)`；[gf.jl:78](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/gf.jl:78) 传入这两个索引，[gf.jl:15](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/gf.jl:15) 将两个变量按这个顺序相乘并除以 $Z_h$。因此取显式末索引 `i=lattice.k=N+1`，正好是式 (58) 的归一化结果。

注意默认返回整段向量的 `occupation(lattice,...)` 在 [occupation.jl:11](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/occupation.jl:11) 遍历 `1:lattice.N`，没有包含 `lattice.k`。此外，在较长轮廓上插入较早节点时，未来的非幺正近似传播不严格相消。因此不能把那种有限步长的中间插入直接等同于较短 $P_h$ 的末时刻占据。本文的单时刻主定理固定末节点；要逐时刻使用这一严格表述，就取相应长度轮廓的末节点。

### 8.3 主定理及等号、近似号的完整分界

在第 2 节有限浴及初态假设下，选 `RealGrassmannLattice1Order` 和预定义 SIAM 精确系统核，所有源单元均为整宽 $h$。若初态系数、IF 系数、MPS 乘积和积分均无额外近似，则

$$
\boxed{
\begin{aligned}
N_O(t)&=\operatorname{Tr}[OW(t)\rho_0W(t)^\dagger],\\
N_{O,h}&=\operatorname{Tr}[OP_h\rho_0P_h^\dagger],
\qquad |N_{O,h}-N_O(t)|\leq C_{t,O}h,\\
N_{O,h}&=\int d\nu_h\,K_{O,h}\mathcal I_h,
\qquad Z_h=\int d\nu_h\,K_h\mathcal I_h=1+O(h),\\
\frac{N_{O,h}}{Z_h}&=N_O(t)+O(h).
\end{aligned}}
\tag{59}
$$

其中传播子近似链已经由式 (26)—(30) 证明；第二行到第三行是式 (35)—(52) 的精确恒等变换；第三行到代码目标由式 (2)—(6)、(55)—(58) 核对。占据数采用实际末节点 `Gt` 插入时，就是 package 已有的观测量构造。一般 $O$ 则须按式 (39) 构造其边界核，不能声称任意一个已有 observable 函数都自动对应它。

因此精确的是近似算符的共同 Grassmann 表示；真实 Hamiltonian 演化与它之间有式 (30)—(34) 的误差，二者不混同。

## 9. 另证 order=1 一般不是二阶：整宽端点的独立反例

只证明 $O(h)$ 上界不排除某些特殊问题更高阶。这里从所有节点整宽的乘积重新计算一个非零一阶主项。

取一个耦合杂质模和一个浴模，其他 SIAM 自旋可置空并解耦，且

$$
A=H_B=0,\qquad B=v(d^\dagger b+b^\dagger d),\quad
\rho_0=|1_d,0_b\rangle\langle1_d,0_b|,\quad O=n_d,\quad v>0.
\tag{60}
$$

该初态是代码能精确构造的单项式 Fock 初态。真实占据为 $\cos^2(vt)$。在单激发子空间，按 $|1,0\rangle,|0,1\rangle$ 排序，单模源指数只有零、一次、二次源项；直接用式 (13) 积分给出精确矩阵

$$
R_h=\begin{pmatrix}1&-ivh\\-ivh&1-v^2h^2/2\end{pmatrix},
\qquad P_h=R_h^{N+1}.
\tag{61}
$$

这里明确是 $N+1$ 次同宽矩阵，而不是改变过端点宽度的另一个乘积。令 $M_0=-iv\sigma_x$、$G=v^2\operatorname{diag}(1,0)$，由矩阵对数展开

$$
\log R_h=hM_0+\frac{h^2}{2}G+O(h^3),
$$

$$
(N+1)\log R_h=tM_0+hM_0+\frac{th}{2}G+O(h^2).
\tag{62}
$$

额外的 $hM_0$ 正是多出一个整宽端点的效应。记 $c=\cos(vt)$、$s=\sin(vt)$，对式 (62) 的矩阵指数作一阶展开，作用于初态得到两个振幅

$$
x=e^{hv^2t/4}\left(c-\frac{3hv}{4}s\right)+O(h^2),\qquad
y=-ie^{hv^2t/4}(s+hvc)+O(h^2).
\tag{63}
$$

例如这可由 $M_0+(h/2)G=(hv^2/4)I-iv\sigma_x+(hv^2/4)\sigma_z$ 的二阶矩阵指数直接展开验证；末尾额外的 $hM_0$ 再分别贡献 $-hvs$ 和 $hvc$。

除以算法自身的分母 $|x|^2+|y|^2$ 后，

$$
\boxed{
n_{d,h}(t)-\cos^2(vt)
=h\left[-2v\cos(vt)\sin(vt)
+\frac v2\cos(vt)\sin^3(vt)\right]+O(h^2).}
\tag{64}
$$

其中第一项来自额外端点时间，第二项来自源重建的非幺正差额。对一般时间该系数非零；例如 $v=0.7,t=1$ 时约为 $-0.618243691187$。所以当前 order=1 离散目标即使对归一化物理占据也一般不是二阶。这个系数由式 (61)—(63) 独立推出，没有借用其他端点布局的主误差。

## 10. 通用 ImpurityHamiltonian 接口另行说明

[generalimpurity/sysdynamics.jl:88](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/generalimpurity/sysdynamics.jl:88) 对各个 Hamiltonian 项应用短时 Grassmann 指数，再由 [bulkconnection.jl:24](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/gvconnections/bulkconnection.jl:24) 添加单位重叠连接。记正常序 Hamiltonian 的符号为 $A_{\rm n}(\bar x,y)$，其合成核是

$$
k_{h,\rm gen}^+(\bar x,y)
=e^{\bar x y}\{1-ihA_{\rm n}(\bar x,y)+O(h^2)\}
=\langle\bar x|(I-ihA)|y\rangle+O(h^2).
\tag{65}
$$

各项对应的 Grassmann 多项式均为偶，所以逐项指数相乘等于符号之和的指数；这不等于相应量子算符指数的乘积。相干态核到有限维算符的线性映射可逆，故唯一对应某个杂质算符 $U_{h,\rm gen}^+=I-ihA+O(h^2)$。例如单模能级的链接系数为 $1-ih\epsilon$，不是 $e^{-ih\epsilon}$。

反分支由代码单独产生，记其对应算符为 $U_{h,\rm gen}^-=I+ihA+O(h^2)$。对标准自伴、正常序条目及伴随相容的实系数实现，它是正分支的伴随；若自定义复杂条目，必须检查输入条目及代码反分支是否保持这个性质，不能仅凭“短时核”就默认成立。

若伴随相容成立，将全部 $U_h$ 换为 $U_{h,\rm gen}^+$，第 5—7 节的精确表示仍成立。单步系统误差为 $O(h^2)$，伸缩求和再增加全局 $O(h)$，所以收敛结论保留，但第 9 节以外的一般主误差要计入这个额外系统误差。

即便只知道两个方向各自的一阶一致性，也可分别建立正、反近似乘积对 $e^{-itH}$、$e^{+itH}$ 的范数界，从而控制其夹心结果的误差；此时不能把它不经检查写成同一个 $P_h\rho_0P_h^\dagger$，也不能直接援引其正性。为保持主定理中每个等号的含义，本文始终选用已核对的预定义 SIAM 核。

## 11. 有限浴结论、连续谱输入和数值审计的边界

有限浴部分已经独立闭合：辅助重建代数、局域差额、包含额外端点的全局误差、完整轮廓拼接、全阶 Gaussian 重求和，以及每个 QuAPI 系数区域都有明确推导。

QuAPI 实际对给定谱函数求积；有限能级求和与谱积分使用相同的矩形、三角权重。这一系数对应性并不自动证明无限浴的动力学收敛定理。若要把式 (59) 的算符范数界推广到连续浴，必须另外控制浴离散极限、耦合算符的有界性及所用常数；本文不把这一步当成已完成的结论。谱求积容差、初态阈值、MPS 截断和记忆截断也都不属于式 (59) 的时间误差证明。

工作区新增审计脚本 [verify_order1_reconstruction.py](/home/zjsun/Desktop/Agent/try_sth/GTEMPO_theory/verify_order1_reconstruction.py)，只使用整宽 order=1 布局，检查如下：

1. 从每个有限维源指数直接提取 Grassmann 系数，保持浴矩阵和外部源的分次乘法，计算完整的 $E_N^+\cdots E_0^+\rho_BE_0^-\cdots E_N^-$ 浴迹；与式 (52) 的全部系数比较，包含最高十二次源项。
2. 比较单模和含相互作用双模的源重建算符演化与 $K\mathcal I$ 积分。
3. 核对末节点的代码式 $a_N^+\bar a_N^-$ 占据插入与辅助迹节点插入一致。
4. 检查式 (64) 的归一化占据误差除以 $h$ 收敛到这里重新计算的整宽端点系数。

结果：12 组完整 IF 多项式比较的最大系数差约为 $2.24\times10^{-16}$；27 组单模与 6 组相互作用双模积分比较的最大差约为 $1.11\times10^{-15}$；式 (64) 的最细步长结果为 $-0.618195360581$，趋近预测的 $-0.618243691187$。

这些只是代数审计，未运行 Julia package，不能代替第 3—9 节的一般证明。源码位置按本次读取的本地版本给出；论文中的阶数说法没有被用作任何关键等号的依据。
