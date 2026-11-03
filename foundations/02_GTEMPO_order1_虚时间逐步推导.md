# 虚时间 GTEMPO 的逐步推导

本文从有限杂质与有限自由 bath 的真实平衡态目标出发，逐节建立离散化及其 Grassmann 积分表示。各节区分三类步骤：辅助对象的精确定义、代数恒等变换，以及带误差估计的实际近似。正文保留推导所需的定义和等式；较长的误差估计、相干态代数与符号依据放在文末附录。

本文的主线是：真实平衡态目标 → 有误差控制的近似传播子 → 精确的共同 Grassmann 积分 → 精确的单元积分 IF。各节先给出结论与证明路线，再展开计算。最后核对当前虚时间 package；代码中的小能量分支差异另行指出，不将其包装为数学恒等式。主定理先使用精确系统核，通用短时系统核在第六节另行处理。

## 1. 有限杂质问题、平衡态目标与一阶虚时间离散化

**本节得到完整虚时间传播的一阶分裂近似 $U_{\mathrm{split}}(\beta)=[U_S(\delta\tau)U_B(\delta\tau)]^N$，并说明由它计算配分函数和观测量的误差为 $O(\delta\tau)$。** 先定义模型和真实平衡态目标，再固定时间步与乘积顺序；虚时间传播的稳定性和误差证明集中在附录 B。本节只做 Lie 分裂，两个单步传播子本身均精确定义。

本节的关系是

$$
\boxed{
\langle O\rangle_\beta
=\frac{\operatorname{Tr}\!\left[O(U_SU_B)^N\right]}
{\operatorname{Tr}\!\left[(U_SU_B)^N\right]}+O(\delta\tau),
\qquad \delta\tau=\beta/N.
}
$$

这里 $U_S=e^{-\delta\tau H_S}$、$U_B=e^{-\delta\tau H_B}$；下面先给出 $H_S,H_B$ 和物理目标的定义，再证明分裂及误差的含义。

### 1.1 模型与真实物理目标

取自然单位 $\hbar=1$，将总 Hamiltonian 写成

$$
\boxed{
H=H_S+H_B,\qquad
H_B=H_{\mathrm{bath}}+H_{\mathrm{hyb}}.
}
\tag{1}
$$

$H_S$ 只作用于杂质，可以包含杂质内部的相互作用；$H_{\mathrm{bath}}$ 描述自由 bath；$H_{\mathrm{hyb}}$ 描述杂质与 bath 的杂化。因此，$H_B$ 是杂质与 bath 的联合算符。

设杂质有 $r$ 个费米模，$\alpha=1,\ldots,r$。用 $d_\alpha$ 表示其湮灭算符，用 $b_k$ 表示有限个 bath 模的湮灭算符，定义

$$
\begin{aligned}
H_{\mathrm{bath}}&=\sum_k\epsilon_k b_k^\dagger b_k,
\qquad F_\alpha=\sum_k V_{\alpha k}b_k,\\
H_{\mathrm{hyb}}&=\sum_\alpha
\left(d_\alpha^\dagger F_\alpha+F_\alpha^\dagger d_\alpha\right).
\end{aligned}
\tag{2}
$$

本文主线采用**各杂质模耦合彼此独立的 bath**：$k$ 统一枚举所有 bath 模，不同 $\alpha$ 的非零 $V_{\alpha k}$ 支持互不相交。因此跨杂质模的 bath 收缩为零；$H_S$ 中仍可有模间相互作用。保留 $\alpha,\gamma$ 两个指标，是为了显式显示这一对角结构。

全部物理模共同满足费米正则反对易关系，包括 impurity mode 与 bath mode 之间的反对易关系。这里的 bath 不含配对项并守粒子数，后文称为正常态自由 bath。本文先固定有限个模，故联合 Fock 空间有限维；$H_S$、$H_B$、$H$ 均为自伴、宇称偶算符。后面的局部积分代数将先用单 impurity mode 展示；本节的分裂估计直接适用于多模以及含相互作用的 $H_S$。

给定逆温 $\beta>0$，设 $O$ 是杂质上的宇称偶、自伴观测量，省略它在 bath 空间上的单位算符。真实的配分函数、带观测量插入的迹与平衡态期望为

$$
\boxed{
\begin{aligned}
Z(\beta)&=\operatorname{Tr}_{S,\mathrm{bath}}e^{-\beta H},\\
Z_O(\beta)&=\operatorname{Tr}_{S,\mathrm{bath}}
\left[Oe^{-\beta H}\right],\\
\langle O\rangle_\beta&=\frac{Z_O(\beta)}{Z(\beta)}.
\end{aligned}
}
\tag{3}
$$

因为 $H$ 自伴且空间有限维，$e^{-\beta H}$ 正定，故 $Z(\beta)>0$。这里的热态是完整联合系统的 $\rho_\beta=e^{-\beta H}/Z(\beta)$，并未假定它在杂质与 bath 之间因子化，也没有另行指定初始密度矩阵。

后续若引入自由 bath 的热态，它必须由 bath 的自由传播权重提取，并同时保留相应的自由 bath 配分函数；它不是本节额外施加的初态假设。

### 1.2 用一阶分裂近似完整虚时间传播

令 $\delta\tau=\beta/N$，其中 $N$ 为正整数。分别定义杂质部分和含杂化的 bath 部分的精确单步传播子：

$$
U_S(\delta\tau)=e^{-\delta\tau H_S},\qquad
U_B(\delta\tau)=e^{-\delta\tau H_B}.
\tag{4}
$$

$U_B$ 作用于杂质与 bath 的联合空间。选择每步先作用 $U_B$、再作用 $U_S$，一阶 Lie 分裂给出

$$
e^{-\delta\tau H}=U_S(\delta\tau)U_B(\delta\tau)+O((\delta\tau)^2).
\tag{5}
$$

这里的 $O((\delta\tau)^2)$ 表示单步算符范数误差。将 $N$ 步相乘，得到**分裂近似传播子**：

$$
\boxed{
U_{\mathrm{split}}(\beta)
=[U_S(\delta\tau)U_B(\delta\tau)]^N,
\qquad N\delta\tau=\beta.
}
\tag{6}
$$

下标 `split` 说明它只做了系统与耦合部分的分裂，尚未引入 Grassmann 积分。全文中这类完整传播子都以 $\beta$ 为参数，单步传播子则以 $\delta\tau$ 为参数；它们均依赖所选步长。对固定有限模型和固定 $\beta$，

$$
\|e^{-\beta H}-U_{\mathrm{split}}(\beta)\|=O(\delta\tau).
\tag{7}
$$

因此，真实配分函数与观测量可写为

$$
\begin{aligned}
Z(\beta)&=\operatorname{Tr}U_{\mathrm{split}}(\beta)+O(\delta\tau),\\
Z_O(\beta)&=\operatorname{Tr}[O U_{\mathrm{split}}(\beta)]+O(\delta\tau),\\
\langle O\rangle_\beta
&=\frac{\operatorname{Tr}[O U_{\mathrm{split}}(\beta)]}
{\operatorname{Tr}U_{\mathrm{split}}(\beta)}+O(\delta\tau).
\end{aligned}
\tag{8}
$$

未标空间的 $\operatorname{Tr}$ 均为杂质与 bath 的联合迹。式 (7)、(8) 是 $N\to\infty$、$N\delta\tau=\beta$ 固定时的上界；它们的虚时间稳定性和归一化证明见附录 B。

## 2. 用 Grassmann 积分近似耦合部分的传播子

**本节要证明：可以将每个 $U_B(\delta\tau)$ 换成一个 Grassmann 积分所定义的近似传播子 $\widetilde U_B(\delta\tau)$，单步差额为 $O((\delta\tau)^2)$，完整传播与观测量的误差仍为 $O(\delta\tau)$。**

本节连接前后两步的核心关系为

$$
\boxed{
U_B(\delta\tau)=\widetilde U_B(\delta\tau)+O((\delta\tau)^2)
\quad\Longrightarrow\quad
[U_S\widetilde U_B]^N=e^{-\beta H}+O(\delta\tau).
}
$$

证明分两步：先把杂化 $H_{\mathrm{hyb}}$ 中的杂质算符换成 Grassmann 变量，并把所得 bath 传播子放入杂质相干态积分；再比较这个积分与真实 $U_B$ 的常数项、一次项。这样的选择使同一对变量同时出现在 bath 传播和杂质相干态中，第三节的共同 $K\times\mathrm{IF}$ 积分就从这里产生。

### 2.1 定义耦合部分的积分近似

先看单杂质模。用独立 Grassmann 变量 $a,\bar a$ 代替 $d,d^\dagger$，定义

$$
\begin{aligned}
H_B^{\mathrm G}(\bar a,a)&=H_{\mathrm{bath}}+\bar a F+F^\dagger a,\\
U_B^{\mathrm G}(\bar a,a;\delta\tau)&=
 e^{-\delta\tau H_B^{\mathrm G}(\bar a,a)}.
\end{aligned}
\tag{9}
$$

上标 $\mathrm G$ 只说明“杂质算符已换成 Grassmann 变量”。$H_B^{\mathrm G}$ 和 $U_B^{\mathrm G}$ 不含杂质算符，仍含 bath 算符；横线是变量标签，$a,\bar a$ 并不是一对复数。它们与全部物理奇算符反对易。这样的外部变量也称为 Grassmann 源，本文不再为“源”引入另一套记号。

现在定义**耦合部分的近似传播子**：

$$
\boxed{
\widetilde U_B(\delta\tau)
=\int D[\bar a,a]\,|a\rangle
U_B^{\mathrm G}(\bar a,a;\delta\tau)\langle\bar a|.
}
\tag{10}
$$

波浪号表示它是对真实 $U_B$ 的近似。与 $U_B^{\mathrm G}$ 不同，$\widetilde U_B$ 是积分后作用于杂质与 bath 的物理算符。此处单模测度为 $D[\bar a,a]=d\bar a\,da\,e^{-\bar a a}$，积分约定为 $\int d\bar a\,da\,a\bar a=1$；多模测度是各模测度之积。相干态约定和逐项积分依据见附录 A。对 $r$ 个杂质模，使用

$$
H_B^{\mathrm G}(\bar a,a)=H_{\mathrm{bath}}
+\sum_\alpha(\bar a_\alpha F_\alpha+F_\alpha^\dagger a_\alpha)
\tag{11}
$$

及同一多模相干态测度。$H_S$ 可含相互作用，它不参与这里的替换。

### 2.2 为什么单步差额只有二阶

为简写式 (10) 这一个积分运算，在附录 A 中记
$Q[f]=\int D[\bar a,a]|a\rangle f\langle\bar a|$。其逐项计算给出 $Q[1]=I$、$Q[H_B^{\mathrm G}]=H_B$，因此

$$
\boxed{
\begin{aligned}
\widetilde U_B(\delta\tau)
&=I-\delta\tau H_B+\frac{(\delta\tau)^2}{2}
Q[(H_B^{\mathrm G})^2]+O((\delta\tau)^3)\\
&=U_B(\delta\tau)+O((\delta\tau)^2).
\end{aligned}
}
\tag{12}
$$

$Q$ 只是分析这个积分的辅助线性映射，不是 package 显式执行的步骤，也不保乘法。式 (10) 是我们的选择，并非 Hamiltonian 唯一规定的选择；选择的依据是式 (12) 的局部一致性和下一节需要的共同变量结构。式 (12) 是对精确积分定义的 Taylor 分析，未把 $\widetilde U_B$ 定义成截断的一次式。

将每步实际替换后，得到后文要表示的**GTEMPO 近似传播子**：

$$
\boxed{
\begin{aligned}
U_{\mathrm{GTEMPO}}(\beta)
&=[U_S(\delta\tau)\widetilde U_B(\delta\tau)]^N,\\
Z_{\mathrm{GTEMPO}}&=\operatorname{Tr}U_{\mathrm{GTEMPO}}(\beta),\\
Z_{O,\mathrm{GTEMPO}}&=\operatorname{Tr}[O U_{\mathrm{GTEMPO}}(\beta)].
\end{aligned}
}
\tag{13}
$$

这几个下标均说明其物理角色：$U_{\mathrm{GTEMPO}}$ 是完整近似传播子，$Z_{\mathrm{GTEMPO}}$ 是近似配分函数，$Z_{O,\mathrm{GTEMPO}}$ 是观测量的近似分子。它们含 $N$ 次系统传播和 $N$ 次耦合传播，总虚时间仍为 $\beta=N\delta\tau$。

固定有限模型时，局部误差可稳定累积为

$$
\boxed{
\begin{aligned}
\|U_{\mathrm{GTEMPO}}(\beta)-e^{-\beta H}\|&=O(\delta\tau),\\
Z_{\mathrm{GTEMPO}}-Z&=O(\delta\tau),\qquad
Z_{O,\mathrm{GTEMPO}}-Z_O=O(\delta\tau),\\
\frac{Z_{O,\mathrm{GTEMPO}}}{Z_{\mathrm{GTEMPO}}}
-\langle O\rangle_\beta&=O(\delta\tau).
\end{aligned}
}
\tag{14}
$$

附录 C 给出二次差额、乘积稳定性及归一化的完整证明。由此，后文可以专注于将 $\operatorname{Tr}[O U_{\mathrm{GTEMPO}}(\beta)]$ 精确写成积分，而不再把近似传播子与真实传播子混同。代码闭合路径的起点与此处乘积顺序的关系留到第六节核对。

## 3. 从算符乘积得到共同变量的系统因子与 bath 权重

**本节直接对 $N$ 步乘积代入积分定义，证明**

$$
\operatorname{Tr}[O(U_S\widetilde U_B)^N]
=\int D[\bar a,a,\bar\chi,\chi]\,K_O\mathcal V_{\mathrm{bath}}.
$$

推导依次为：代入每个 $\widetilde U_B$ → 保持 bath 次序并收集系统核 → 将杂质迹写成相干态积分。系统因子 $K_O$ 与 bath 权重 $\mathcal V_{\mathrm{bath}}$ 由同一次积分自然产生；末端 bra 的负号同时给出反周期闭合。本节全为恒等变换，尚未使用 bath 的 Gaussian 性。

### 3.1 代入全部耦合传播子并收集系统核

以下省略固定的单步参数，$U_S=e^{-\delta\tau H_S}$。为第 $j$ 个耦合单元引入独立变量 $a_j,\bar a_j$，$j=0,\ldots,N-1$；每个变量都包含 $r$ 个杂质模分量。式 (10) 成为

$$
U_{B,j}^{\mathrm G}=e^{-\delta\tau H_B^{\mathrm G}(\bar a_j,a_j)},
\qquad
\widetilde U_B=\int D[\bar a_j,a_j]\,|a_j\rangle
U_{B,j}^{\mathrm G}\langle\bar a_j|.
\tag{15}
$$

所有节点变量彼此反对易，且均与物理奇算符反对易。规范次序先按节点递增、再按模次序排列；稍后出现的迹变量放在这些节点变量之前。每对测度为偶，交换整对测度不产生符号；单个生成元的交换仍须计入符号。

把式 (15) 逐个代入式 (13)，从右向左的节点依次为 $0,1,\ldots,N-1$：

$$
\begin{aligned}
U_{\mathrm{GTEMPO}}(\beta)
&=\int\prod_{j=0}^{N-1}D[\bar a_j,a_j]\,
\big(U_S|a_{N-1}\rangle U_{B,N-1}^{\mathrm G}\langle\bar a_{N-1}|\big)
\cdots
\big(U_S|a_0\rangle U_{B,0}^{\mathrm G}\langle\bar a_0|\big)\\
&=\int\prod_{j=0}^{N-1}D[\bar a_j,a_j]\,
U_S|a_{N-1}\rangle
\left[\prod_{j=1}^{N-1}\langle\bar a_j|U_S|a_{j-1}\rangle\right]
\langle\bar a_0|\,
U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}.
\end{aligned}
\tag{16}
$$

第一行只是积分定义的代入；第二行将 bath 因子向右收集。含 Grassmann 系数的相干 ket、bra 以及 $U_{B,j}^{\mathrm G}$ 作为整体均为偶，且系统与 bath 的物理因子不同，因此允许这种重排，详见附录 D。**不同 bath 因子之间未作交换**，始终保持 $U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}$ 的次序。各系统核为偶 Grassmann 标量，彼此可交换。$N=1$ 时内部乘积按空乘积取 $1$，无需另作推导。

### 3.2 取迹，得到共同积分与反周期闭合

对式 (16) 左乘 $O$ 并取联合迹，直接将杂质迹写成相干态积分：

$$
\begin{aligned}
Z_{O,\mathrm{GTEMPO}}
&=\int\prod_{j=0}^{N-1}D[\bar a_j,a_j]\,
\operatorname{Tr}_S\!\left[
OU_S|a_{N-1}\rangle
\left(\prod_{j=1}^{N-1}\langle\bar a_j|U_S|a_{j-1}\rangle\right)
\langle\bar a_0|\right]
\operatorname{Tr}_B\!\left[U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}\right]\\
&=\int D[\bar\chi,\chi]\prod_{j=0}^{N-1}D[\bar a_j,a_j]\,
\langle-\bar\chi|OU_S|a_{N-1}\rangle
\left[\prod_{j=1}^{N-1}\langle\bar a_j|U_S|a_{j-1}\rangle\right]
\langle\bar a_0|\chi\rangle\,
\operatorname{Tr}_B\!\left[U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}\right].
\end{aligned}
\tag{17}
$$

这里 $\chi,\bar\chi$ 与全部节点变量独立，不进入 bath 指数。求迹的负号及含 Grassmann 系数时的分次符号在附录 D 证明。现在把同一积分中的因子分别命名：

$$
\boxed{
\begin{aligned}
D[\bar a,a,\bar\chi,\chi]
&=D[\bar\chi,\chi]\prod_{j=0}^{N-1}D[\bar a_j,a_j],\\
K_O&=\langle-\bar\chi|OU_S|a_{N-1}\rangle
\left[\prod_{j=1}^{N-1}\langle\bar a_j|U_S|a_{j-1}\rangle\right]
\langle\bar a_0|\chi\rangle,\\
\mathcal V_{\mathrm{bath}}[\bar a,a]
&=\operatorname{Tr}_B\!\left[U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}\right].
\end{aligned}
}
\tag{18}
$$

因此精确得到

$$
\boxed{
\begin{aligned}
Z_{O,\mathrm{GTEMPO}}
&=\int D[\bar a,a,\bar\chi,\chi]\,K_O\mathcal V_{\mathrm{bath}}[\bar a,a],\\
Z_{\mathrm{GTEMPO}}
&=\int D[\bar a,a,\bar\chi,\chi]\,K_I\mathcal V_{\mathrm{bath}}[\bar a,a].
\end{aligned}
}
\tag{19}
$$

同一节点的 $a_j,\bar a_j$ 既出现在 $U_{B,j}^{\mathrm G}$ 中，又分别进入相邻系统核；这是式 (15) 的同一次积分定义的结果，并未将两套独立变量事后认同。内部有 $N-1$ 个 $U_S$，末端 $OU_S$ 中还有一个，共 $N$ 个；$\chi$ 不增加传播时间，也不代表初态。

末端负号已经包含反周期边界。具体地，先积分迹变量可得

$$
\boxed{
\int D[\bar\chi,\chi]\,
\langle-\bar\chi|OU_S|a_{N-1}\rangle
\langle\bar a_0|\chi\rangle
=\langle-\bar a_0|OU_S|a_{N-1}\rangle.
}
\tag{20}
$$

例如单模 $OU_S=I$ 时，两侧均为 $e^{-\bar a_0a_{N-1}}$。这就是费米路径绕虚时间圈闭合时接缝取负的离散表达；最早和最晚节点仍是独立积分变量。相应的逐项符号核验见附录 D。

## 4. 离散 Grassmann 变量如何进入连续时间有序积分

**本节证明：有限个 bath 指数的乘积，可以精确写成一个有序指数；其中仍然只有原来的 $N$ 组离散 Grassmann 变量。** 结论的完整写法是

$$
\boxed{
U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}
=e^{-\beta H_{\mathrm{bath}}}\mathcal T\exp\!\left[
-\sum_{j=0}^{N-1}\int_{C_j}d\tau\,
\sum_\alpha\big(\bar a_{j\alpha}F_\alpha(\tau)
+F_\alpha^\dagger(\tau)a_{j\alpha}\big)\right].
}
$$

这里 $C_j=[j\delta\tau,(j+1)\delta\tau)$。右侧每个 $a_{j\alpha},\bar a_{j\alpha}$ 都显式保留。连续时间 $\tau$ 是展开有限指数时出现的积分变量；它既不是新的 Grassmann 自由度，也不要求取 $N\to\infty$。

证明分三步：从一个单元的指数推出精确 Dyson 展开 → 用两个单元说明时间移位和离散变量的保留 → 拼接 $N$ 个单元，并从自由传播中提取 $Z_{\mathrm{bath}}$。本节没有新增时间近似。

### 4.1 单个指数的精确展开：连续积分变量从哪里来

对第 $j$ 个单元，式 (15) 中的指数是

$$
\begin{aligned}
H_{\mathrm{hyb},j}^{\mathrm G}
&=\sum_\alpha\big(\bar a_{j\alpha}F_\alpha+F_\alpha^\dagger a_{j\alpha}\big),\\
U_{B,j}^{\mathrm G}
&=e^{-\delta\tau(H_{\mathrm{bath}}+H_{\mathrm{hyb},j}^{\mathrm G})}.
\end{aligned}
\tag{21}
$$

单元内 $a_{j\alpha},\bar a_{j\alpha}$ 固定不变。为了分出自由 bath 传播，定义

$$
\begin{aligned}
U_{I,j}(u)
&=e^{uH_{\mathrm{bath}}}e^{-u(H_{\mathrm{bath}}+H_{\mathrm{hyb},j}^{\mathrm G})},\\
H_{\mathrm{hyb},j}^{\mathrm G}(u)
&=e^{uH_{\mathrm{bath}}}H_{\mathrm{hyb},j}^{\mathrm G}e^{-uH_{\mathrm{bath}}},
\qquad 0\leq u\leq\delta\tau.
\end{aligned}
\tag{22}
$$

这是有限指数的相互作用绘景分解。直接求导再积分，得到

$$
\begin{aligned}
U_{I,j}'(u)&=-H_{\mathrm{hyb},j}^{\mathrm G}(u)U_{I,j}(u),
\qquad U_{I,j}(0)=I,\\
U_{I,j}(u)&=I-\int_0^u dv\,H_{\mathrm{hyb},j}^{\mathrm G}(v)U_{I,j}(v).
\end{aligned}
\tag{23}
$$

向右侧反复代入这个积分方程，就得到

$$
\begin{aligned}
U_{I,j}(\delta\tau)={}&I
-\int_0^{\delta\tau}du_1\,H_{\mathrm{hyb},j}^{\mathrm G}(u_1)\\
&+\int_0^{\delta\tau}du_1\int_0^{u_1}du_2\,
H_{\mathrm{hyb},j}^{\mathrm G}(u_1)H_{\mathrm{hyb},j}^{\mathrm G}(u_2)\\
&-\int_0^{\delta\tau}du_1\int_0^{u_1}du_2\int_0^{u_2}du_3\,
H_{\mathrm{hyb},j}^{\mathrm G}(u_1)H_{\mathrm{hyb},j}^{\mathrm G}(u_2)
H_{\mathrm{hyb},j}^{\mathrm G}(u_3)+\cdots .
\end{aligned}
\tag{24}
$$

这就是 Dyson 展开。$u_1,u_2,\ldots$ 来自求解式 (23) 的积分方程，不来自额外的时间离散化。每个因子含一个本单元的 Grassmann 生成元；一个单元只有 $2r$ 个生成元，故超过 $2r$ 个因子的项为零。上式的省略号包含全部剩余非零项，不表示截断误差。

相似变换把 $F_\alpha$ 变为 $F_\alpha(u)=e^{uH_{\mathrm{bath}}}F_\alpha e^{-uH_{\mathrm{bath}}}$，把 $F_\alpha^\dagger$ 也按同样方式变换；具体自由模表达将在下一小节写出。用 $\mathcal T$ 表示较大时间在左，将完整展开记为有序指数，便有

$$
\boxed{
U_{B,j}^{\mathrm G}
=e^{-\delta\tau H_{\mathrm{bath}}}\mathcal T\exp\!\left[
-\int_0^{\delta\tau}du\,
\sum_\alpha\big(\bar a_{j\alpha}F_\alpha(u)
+F_\alpha^\dagger(u)a_{j\alpha}\big)\right].
}
\tag{25}
$$

$\mathcal T$ 排序整个宇称偶因子 $H_{\mathrm{hyb},j}^{\mathrm G}(u)$，交换这类因子的时间标签时不额外附加费米符号。单个 Grassmann 变量和 bath 奇算符的交换仍须按分次规则处理。式 (25) 是初值问题的精确解，不能去掉 $\mathcal T$ 后把它当作普通指数的拆分。

### 4.2 两个单元：积分时间变化，离散变量仍留在各自的单元

相似变换只改变 bath 算符的时间标签，Grassmann 变量不变。具体地，

$$
\begin{aligned}
F_\alpha(\tau)
&=e^{\tau H_{\mathrm{bath}}}F_\alpha e^{-\tau H_{\mathrm{bath}}}
=\sum_kV_{\alpha k}e^{-\epsilon_k\tau}b_k,\\
F_\alpha^\dagger(\tau)
&=e^{\tau H_{\mathrm{bath}}}F_\alpha^\dagger e^{-\tau H_{\mathrm{bath}}}
=\sum_kV_{\alpha k}^*e^{+\epsilon_k\tau}b_k^\dagger.
\end{aligned}
\tag{26}
$$

**$F_\alpha^\dagger(\tau)$ 表示产生算符的绘景像，一般不是 $F_\alpha(\tau)$ 的厄米伴随。** 后者的指数仍为 $e^{-\epsilon_k\tau}$。后续始终使用上式的相似变换定义。

对两个单元，先不合并有序指数，逐步写为

$$
\begin{aligned}
U_{B,1}^{\mathrm G}U_{B,0}^{\mathrm G}
&=e^{-\delta\tau H_{\mathrm{bath}}}U_{I,1}(\delta\tau)
 e^{-\delta\tau H_{\mathrm{bath}}}U_{I,0}(\delta\tau)\\
&=e^{-2\delta\tau H_{\mathrm{bath}}}
\big[e^{\delta\tau H_{\mathrm{bath}}}U_{I,1}(\delta\tau)
 e^{-\delta\tau H_{\mathrm{bath}}}\big]U_{I,0}(\delta\tau).
\end{aligned}
\tag{27}
$$

对式 (24) 中 $j=1$ 的每一项作方括号里的相似变换。自由 Hamiltonian 为偶，与外部变量交换，所以每个因子变为

$$
e^{\delta\tau H_{\mathrm{bath}}}H_{\mathrm{hyb},1}^{\mathrm G}(u)
 e^{-\delta\tau H_{\mathrm{bath}}}
=\sum_\alpha\big[\bar a_{1\alpha}F_\alpha(\delta\tau+u)
+F_\alpha^\dagger(\delta\tau+u)a_{1\alpha}\big].
$$

令 $\tau=\delta\tau+u$，两个单元的乘积便是

$$
\boxed{
\begin{aligned}
U_{B,1}^{\mathrm G}U_{B,0}^{\mathrm G}
={}&e^{-2\delta\tau H_{\mathrm{bath}}}\\
&\times\mathcal T\exp\!\left[-\int_{\delta\tau}^{2\delta\tau}d\tau\,
\sum_\alpha\big(\bar a_{1\alpha}F_\alpha(\tau)
+F_\alpha^\dagger(\tau)a_{1\alpha}\big)\right]\\
&\times\mathcal T\exp\!\left[-\int_0^{\delta\tau}d\tau\,
\sum_\alpha\big(\bar a_{0\alpha}F_\alpha(\tau)
+F_\alpha^\dagger(\tau)a_{0\alpha}\big)\right].
\end{aligned}
}
\tag{28}
$$

这时可以直接看到：第二单元的指数仍含 $a_1,\bar a_1$，第一单元的指数仍含 $a_0,\bar a_0$。只改变了展开 bath 算符所用的连续时间标签，没有积分掉或合并任何 Grassmann 变量。

为了将这种表达写成整个区间上的一个积分，定义**原有离散变量的分段常表示**：

$$
\boxed{
\begin{aligned}
C_j&=[j\delta\tau,(j+1)\delta\tau),\\
a_\alpha(\tau)&=a_{j\alpha},\qquad
\bar a_\alpha(\tau)=\bar a_{j\alpha},\qquad \tau\in C_j,\\
H_{\mathrm{hyb}}^{\mathrm G}(\tau)
&=\sum_\alpha\big[\bar a_\alpha(\tau)F_\alpha(\tau)
+F_\alpha^\dagger(\tau)a_\alpha(\tau)\big].
\end{aligned}
}
\tag{29}
$$

例如只考虑前两个单元时，$a(\tau)$ 在 $[0,\delta\tau)$ 上就是 $a_0$，在 $[\delta\tau,2\delta\tau)$ 上就是 $a_1$。这不是给每个连续时间点引入一对新变量，而是给同一组离散变量添加单元标签。$H_{\mathrm{hyb}}^{\mathrm G}(\tau)$ 因此同时含连续时间的 bath 绘景算符与分段常的离散源。

现在检验两个有序指数如何合并。将式 (28) 的两个指数各自展开，乘积的零、一次、二次项依次为

$$
\begin{aligned}
I&-\left[\int_0^{\delta\tau}d\tau
+\int_{\delta\tau}^{2\delta\tau}d\tau\right]
H_{\mathrm{hyb}}^{\mathrm G}(\tau)\\
&+\int_0^{\delta\tau}d\tau\int_0^\tau d\tau'\,
H_{\mathrm{hyb}}^{\mathrm G}(\tau)H_{\mathrm{hyb}}^{\mathrm G}(\tau')\\
&+\int_{\delta\tau}^{2\delta\tau}d\tau\int_{\delta\tau}^\tau d\tau'\,
H_{\mathrm{hyb}}^{\mathrm G}(\tau)H_{\mathrm{hyb}}^{\mathrm G}(\tau')\\
&+\int_{\delta\tau}^{2\delta\tau}d\tau\int_0^{\delta\tau}d\tau'\,
H_{\mathrm{hyb}}^{\mathrm G}(\tau)H_{\mathrm{hyb}}^{\mathrm G}(\tau')+\cdots .
\end{aligned}
\tag{30}
$$

二阶的两块三角形来自单元内部，矩形来自第二单元一次项乘第一单元一次项。矩形中的乘积显式包含

$$
\sum_{\alpha\gamma}
\big[\bar a_{1\alpha}F_\alpha(\tau)+F_\alpha^\dagger(\tau)a_{1\alpha}\big]
\big[\bar a_{0\gamma}F_\gamma(\tau')+F_\gamma^\dagger(\tau')a_{0\gamma}\big].
$$

三块恰好覆盖 $0\leq\tau'\leq\tau\leq2\delta\tau$；原有两组变量在每块中各取其单元的值。因此乘积等于整个两单元区间的 Dyson 展开。

### 4.3 拼接所有单元，并从自由传播中提取 bath 配分函数

全阶拼接使用同一个规则：把全区间中一组有序时间，按落入各单元的个数分类。每一类恰好对应各单元 Dyson 项的一个乘积；较晚单元整体在左，单元内部再按时间排序。分类不重不漏，单元边界为零测集。这对有限 $N$ 的每一阶都成立。

因而精确得到

$$
\boxed{
\begin{aligned}
U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}
&=e^{-\beta H_{\mathrm{bath}}}\mathcal T\exp\!\left[
-\sum_{j=0}^{N-1}\int_{C_j}d\tau\,
\sum_\alpha\big(\bar a_{j\alpha}F_\alpha(\tau)
+F_\alpha^\dagger(\tau)a_{j\alpha}\big)\right]\\
&=e^{-\beta H_{\mathrm{bath}}}\mathcal T\exp\!\left[
-\int_0^\beta d\tau\,H_{\mathrm{hyb}}^{\mathrm G}(\tau)\right].
\end{aligned}
}
\tag{31}
$$

**第一行把所有离散源写出来，第二行只是用式 (29) 的分段常函数简写。** 两行是同一个外代数中的同一个 bath 算符多项式。$\int d\tau$ 只积分连续时间的系数函数；它不执行 Grassmann 积分。真正对这些离散源的积分仍是第三节的 $\int D[\bar a,a,\bar\chi,\chi]$。

从左侧的自由传播权重提取

$$
Z_{\mathrm{bath}}=\operatorname{Tr}_B e^{-\beta H_{\mathrm{bath}}}
=\prod_k(1+e^{-\beta\epsilon_k}),\qquad
\rho_{\mathrm{bath}}^{\mathrm{free}}
=\frac{e^{-\beta H_{\mathrm{bath}}}}{Z_{\mathrm{bath}}}.
\tag{32}
$$

这只是已有权重的归一化，没有假定完整平衡态是因子化初态。定义 $\langle X\rangle_{\mathrm{bath}}=\operatorname{Tr}_B(\rho_{\mathrm{bath}}^{\mathrm{free}}X)$，并在下一节明确它对含源多项式的作用，式 (18) 的 bath 权重便是

$$
\boxed{
\begin{aligned}
\mathcal V_{\mathrm{bath}}[\bar a,a]
&=Z_{\mathrm{bath}}\,\mathrm{IF}[\bar a,a],\\
\mathrm{IF}[\bar a,a]
&=\left\langle\mathcal T\exp\!\left[
-\int_0^\beta d\tau\,H_{\mathrm{hyb}}^{\mathrm G}(\tau)\right]
\right\rangle_{\mathrm{bath}}.
\end{aligned}
}
\tag{33}
$$

自由传播始终在最左，推导没有对含 Grassmann 奇系数的 bath 迹使用普通循环性。$Z_{\mathrm{bath}}$ 会保留在配分函数与观测量分子中，直到取比值时才约去。

## 5. 逐项做 bath 平均，再重求和为 QuAPI 影响泛函

**本节计算式 (33) 的 bath 平均，证明**

$$
\boxed{
\mathrm{IF}=1+\Phi_{\mathrm{bath}}
+\frac{\Phi_{\mathrm{bath}}^2}{2!}+\cdots
=e^{\Phi_{\mathrm{bath}}},\qquad
\Phi_{\mathrm{bath}}
=\sum_{i\alpha,j\gamma}\bar a_{i\alpha}\Delta_{i\alpha,j\gamma}a_{j\gamma}.
}
$$

这里 $\Phi_{\mathrm{bath}}$ 将由实际的二阶计算确定。证明依次回答：bath 平均去掉什么、保留什么 → 为什么只有偶数阶 → 二阶的两个符号是什么 → 四阶三种配对为何给出 $\Phi_{\mathrm{bath}}^2/2!$ → 一般偶数阶如何组成指数 → 如何提出分段常变量，得到 QuAPI 单元积分。本节直到最后与真实物理目标比较之前，均为精确定义和恒等变换。

### 5.1 bath 平均只计算 bath 系数，Grassmann 变量仍然保留

式 (33) 的平均只在有限 bath Fock 空间中取迹。先按既定分次规则把每个源单项式放左、bath 算符系数放右，再逐项定义

$$
\boxed{
\left\langle\sum_m m(\bar a,a)X_m\right\rangle_{\mathrm{bath}}
=\sum_m m(\bar a,a)\,
\operatorname{Tr}_B(\rho_{\mathrm{bath}}^{\mathrm{free}}X_m).
}
\tag{34}
$$

这里 $m$ 遍历源单项式，$X_m$ 是其 bath 算符系数。例如，

$$
\left\langle\bar a_{i\alpha}a_{j\gamma}
F_\alpha(\tau)F_\gamma^\dagger(\tau')\right\rangle_{\mathrm{bath}}
=\bar a_{i\alpha}a_{j\gamma}
\left\langle F_\alpha(\tau)F_\gamma^\dagger(\tau')\right\rangle_{\mathrm{bath}}.
$$

右侧的 bath 平均是复数系数，$\bar a_{i\alpha}a_{j\gamma}$ 原样留下。因此 IF 是原有离散 Grassmann 变量的多项式；它不再含 bath 算符，但尚未完成杂质变量的 Berezin 积分。

将有序指数完全展开，得到

$$
\begin{aligned}
\mathrm{IF}
={}&\sum_{p=0}^{2rN}\frac{(-1)^p}{p!}
\int_0^\beta d\tau_1\cdots\int_0^\beta d\tau_p\,
\left\langle\mathcal T
H_{\mathrm{hyb}}^{\mathrm G}(\tau_1)\cdots
H_{\mathrm{hyb}}^{\mathrm G}(\tau_p)\right\rangle_{\mathrm{bath}}.
\end{aligned}
\tag{35}
$$

$p=0$ 项为 $1$。这里把 Dyson 的有序单纯形改写为整个时间立方体：$\mathcal T$ 负责排序，$1/p!$ 补偿同一有序区域的 $p!$ 次覆盖。**完整因子** $H_{\mathrm{hyb}}^{\mathrm G}$ 的总宇称为偶，所以这个改写没有额外置换符号；因子内部的分次符号仍须保留。

每个因子含一个 bath 奇算符，$p$ 个因子的 bath 系数宇称为 $p$ 模 $2$。自由 bath 热态宇称偶。用 bath 宇称算符作相似变换时，热态不变、奇算符变号，故奇算符的平均等于自身的负值，只能为零。于是所有奇数阶消失：

$$
\left\langle\mathcal T
H_{\mathrm{hyb}}^{\mathrm G}(\tau_1)\cdots
H_{\mathrm{hyb}}^{\mathrm G}(\tau_{2m+1})\right\rangle_{\mathrm{bath}}=0.
$$

以下用 $\mathrm{IF}^{(p)}$ 表示 Grassmann 总次数为 $p$ 的部分。于是 $\mathrm{IF}=1+\mathrm{IF}^{(2)}+\mathrm{IF}^{(4)}+\cdots$。上标表示源的次数，不是时间近似的阶数；现在没有丢弃任何非零项。

### 5.2 实际计算二阶项，固定两个时间区域的符号

把二阶部分命名为 $\Phi_{\mathrm{bath}}$，目的是先算出一个具体的二次多项式，再检查高阶是否由它生成。二阶 Dyson 负号相乘为正，因此

$$
\Phi_{\mathrm{bath}}:=\mathrm{IF}^{(2)}
=\int_0^\beta d\tau\int_0^\tau d\tau'\,
\left\langle H_{\mathrm{hyb}}^{\mathrm G}(\tau)
H_{\mathrm{hyb}}^{\mathrm G}(\tau')\right\rangle_{\mathrm{bath}}.
\tag{36}
$$

先用单杂质模看清符号。对 $\tau>\tau'$，两个因子展开成四项：

$$
\begin{aligned}
&(\bar a(\tau)F(\tau)+F^\dagger(\tau)a(\tau))
 (\bar a(\tau')F(\tau')+F^\dagger(\tau')a(\tau'))\\
={}&(\bar a(\tau)F(\tau))(\bar a(\tau')F(\tau'))
 +(\bar a(\tau)F(\tau))(F^\dagger(\tau')a(\tau'))\\
&+(F^\dagger(\tau)a(\tau))(\bar a(\tau')F(\tau'))
 +(F^\dagger(\tau)a(\tau))(F^\dagger(\tau')a(\tau')).
\end{aligned}
$$

第一、四项的 bath 平均分别含 $\langle FF\rangle_{\mathrm{bath}}$、$\langle F^\dagger F^\dagger\rangle_{\mathrm{bath}}$，在守粒子数的自由热态中均为零。另两项按分次规则重排：

$$
\begin{aligned}
(\bar a(\tau)F(\tau))(F^\dagger(\tau')a(\tau'))
&=\bar a(\tau)a(\tau')F(\tau)F^\dagger(\tau'),\\
(F^\dagger(\tau)a(\tau))(\bar a(\tau')F(\tau'))
&=-\bar a(\tau')a(\tau)F^\dagger(\tau)F(\tau').
\end{aligned}
\tag{37}
$$

第一项中，把 $a(\tau')$ 移过两个 bath 奇算符，符号为正。第二项先把两个源提到左侧，再把 $a(\tau)\bar a(\tau')$ 排成 $\bar a(\tau')a(\tau)$，留下一个负号。这些符号来自代数重排，没有更改 bath 迹的次序。

所需的两个纯 bath 平均可直接由自由模占据 $n_k=(1+e^{\beta\epsilon_k})^{-1}$ 计算。写出多杂质模形式，单模时省略 $\alpha,\gamma$：

$$
\begin{aligned}
C^>_{\alpha\gamma}(\tau,\tau')
&:=\left\langle F_\alpha(\tau)F_\gamma^\dagger(\tau')\right\rangle_{\mathrm{bath}}
=\sum_kV_{\alpha k}V_{\gamma k}^*(1-n_k)e^{-\epsilon_k(\tau-\tau')},\\
C^<_{\alpha\gamma}(\tau,\tau')
&:=\left\langle F_\gamma^\dagger(\tau')F_\alpha(\tau)\right\rangle_{\mathrm{bath}}
=\sum_kV_{\alpha k}V_{\gamma k}^*n_ke^{-\epsilon_k(\tau-\tau')}.
\end{aligned}
\tag{38}
$$

因此对 $\tau>\tau'$，完整二因子平均为

$$
\begin{aligned}
&\left\langle H_{\mathrm{hyb}}^{\mathrm G}(\tau)
H_{\mathrm{hyb}}^{\mathrm G}(\tau')\right\rangle_{\mathrm{bath}}\\
&\quad=\sum_{\alpha\gamma}\Big[
\bar a_\alpha(\tau)C^>_{\alpha\gamma}(\tau,\tau')a_\gamma(\tau')
-\bar a_\alpha(\tau')C^<_{\alpha\gamma}(\tau',\tau)a_\gamma(\tau)\Big].
\end{aligned}
\tag{39}
$$

在本文独立 bath 的假设下，$C^{>,<}_{\alpha\gamma}=0$（$\alpha\ne\gamma$）。这里保留双指标以显示多模推广；$H_S$ 的相互作用始终位于第三节的系统因子中，不影响这一步 bath 平均。

将第二项中的两个积分变量互换名称，两个三角区域就成为整个时间方形。定义连续时间杂化核：

$$
\boxed{
\begin{aligned}
\Delta_{\alpha\gamma}(\tau,\tau')
&=\begin{cases}
C^>_{\alpha\gamma}(\tau,\tau'),&\tau>\tau',\\
-C^<_{\alpha\gamma}(\tau,\tau'),&\tau<\tau',
\end{cases}\\
\Phi_{\mathrm{bath}}
&=\sum_{\alpha\gamma}\int_0^\beta d\tau\int_0^\beta d\tau'\,
\bar a_\alpha(\tau)\Delta_{\alpha\gamma}(\tau,\tau')a_\gamma(\tau').
\end{aligned}
}
\tag{40}
$$

这已经算出了 $\Phi_{\mathrm{bath}}$，但尚未说明 IF 为什么是它的指数。等时线的取值不影响积分；同一离散单元中的两个三角区域则都有有限面积，必须保留。

### 5.3 四阶的三种配对：为什么恰好是二阶项的平方除以二

四阶来自四个杂化因子的 bath 平均。先固定 $\tau_1>\tau_2>\tau_3>\tau_4$，用单模展开中的一个具体乘积检查 Wick 符号。以下仅在本段记 $F_i=F(\tau_i)$、$F_i^\dagger=F^\dagger(\tau_i)$；$i$ 在这里标记连续积分时间。源仍写作 $a(\tau_i),\bar a(\tau_i)$，不另行编号；若两个时间落在同一单元，它们就取同一个原有离散生成元。

把源提出后，这个具体乘积为

$$
\begin{aligned}
&\left\langle(\bar a(\tau_1)F_1)(\bar a(\tau_2)F_2)
(F_3^\dagger a(\tau_3))(F_4^\dagger a(\tau_4))\right\rangle_{\mathrm{bath}}\\
&\quad=\bar a(\tau_1)\bar a(\tau_2)a(\tau_3)a(\tau_4)\,
\left\langle F_1F_2F_3^\dagger F_4^\dagger\right\rangle_{\mathrm{bath}}\\
&\quad=\bar a(\tau_1)\bar a(\tau_2)a(\tau_3)a(\tau_4)
\Big[-C^>(\tau_1,\tau_3)C^>(\tau_2,\tau_4)
+C^>(\tau_1,\tau_4)C^>(\tau_2,\tau_3)\Big].
\end{aligned}
\tag{41}
$$

纯 bath 的四算符 Wick 公式有三种配对：$(12)(34)$、$(13)(24)$、$(14)(23)$，符号依次为正、负、正。在这个特定乘积里，第一种因 $\langle F_1F_2\rangle_{\mathrm{bath}}=0$ 而消失。另两种也可以直接写成两个完整二次项的乘积：

$$
\begin{aligned}
&\big[\bar a(\tau_1)C^>(\tau_1,\tau_3)a(\tau_3)\big]
 \big[\bar a(\tau_2)C^>(\tau_2,\tau_4)a(\tau_4)\big]\\
&\quad+\big[\bar a(\tau_1)C^>(\tau_1,\tau_4)a(\tau_4)\big]
 \big[\bar a(\tau_2)C^>(\tau_2,\tau_3)a(\tau_3)\big].
\end{aligned}
$$

第一行的源次序是 $\bar a(\tau_1)a(\tau_3)\bar a(\tau_2)a(\tau_4)=-\bar a(\tau_1)\bar a(\tau_2)a(\tau_3)a(\tau_4)$；第二行是 $\bar a(\tau_1)a(\tau_4)\bar a(\tau_2)a(\tau_3)=+\bar a(\tau_1)\bar a(\tau_2)a(\tau_3)a(\tau_4)$。所以在**完整二因子平均之积**的写法里，两种配对前都是加号，纯 bath 的负号已由源排序承担。

对四个完整杂化因子作相同处理：将每个因子的两项以及所有杂质模指标都展开，对每一项应用纯 bath Wick 公式，再把源恢复到对应的配对次序。得到

$$
\begin{aligned}
\left\langle\mathcal T H_1H_2H_3H_4\right\rangle_{\mathrm{bath}}
={}&\left\langle\mathcal T H_1H_2\right\rangle_{\mathrm{bath}}
   \left\langle\mathcal T H_3H_4\right\rangle_{\mathrm{bath}}\\
&+\left\langle\mathcal T H_1H_3\right\rangle_{\mathrm{bath}}
   \left\langle\mathcal T H_2H_4\right\rangle_{\mathrm{bath}}\\
&+\left\langle\mathcal T H_1H_4\right\rangle_{\mathrm{bath}}
   \left\langle\mathcal T H_2H_3\right\rangle_{\mathrm{bath}},
\qquad H_i:=H_{\mathrm{hyb}}^{\mathrm G}(\tau_i).
\end{aligned}
\tag{42}
$$

这对任意时间标签成立；$\mathcal T$ 排序完整偶因子。源符号为何对每一种算符类型和配对都成立，见附录 E.2；有限自由 bath 的 Wick 递推证明见附录 E.1。这里没有假定带源的 bath 算符彼此交换，也没有对它们使用普通迹循环性。

接下来计算组合系数。先把二阶有序积分写成全时间方形：

$$
\Phi_{\mathrm{bath}}
=\frac12\int_0^\beta d\tau\int_0^\beta d\tau'\,
\left\langle\mathcal T H_{\mathrm{hyb}}^{\mathrm G}(\tau)
H_{\mathrm{hyb}}^{\mathrm G}(\tau')\right\rangle_{\mathrm{bath}}.
\tag{43}
$$

四阶使用式 (35) 的整个时间立方体，系数为 $1/4!$。式 (42) 的每一种配对，都把四个积分变量分成互不相交的两对。重命名哑变量后，每一种配对的积分都等于上式中二因子积分的平方，即 $(2\Phi_{\mathrm{bath}})^2$。因此

$$
\boxed{
\mathrm{IF}^{(4)}
=\frac{3}{4!}\,(2\Phi_{\mathrm{bath}})^2
=\frac{\Phi_{\mathrm{bath}}^2}{2!}.
}
\tag{44}
$$

这里“三种配对积分相等”的依据是**整个时间立方体中的变量重命名**。在固定 $\tau_1>\tau_2>\tau_3>\tau_4$ 的有序区域里，不能直接声称三种配对各自的积分相等。

### 5.4 推广到全部偶数阶，得到精确二次指数

现在推广刚才的两个步骤：先配对，再计数。对 $2m$ 个完整杂化因子，纯 bath Wick 定理给出所有两两配对；每个配对的 bath 置换符号，均由相同的源置换吸收。所以每一项是 $m$ 个完整二因子平均之积。较长的分次符号证明放在附录 E.2。

配对数为 $(2m)!/(2^m m!)$：先将 $2m$ 个标签排列，再除去每对内部的两种排列和 $m$ 对之间的排列。在全时间立方体上，每种配对积分重命名后均为 $(2\Phi_{\mathrm{bath}})^m$，于是

$$
\boxed{
\mathrm{IF}^{(2m)}
=\frac{1}{(2m)!}\frac{(2m)!}{2^m m!}
 (2\Phi_{\mathrm{bath}})^m
=\frac{\Phi_{\mathrm{bath}}^m}{m!}.
}
\tag{45}
$$

四阶的 $3/4!$ 正是这个公式在 $m=2$ 时的特例。包括零阶以及已经消失的奇数阶，便得到

$$
\boxed{
\mathrm{IF}[\bar a,a]
=\sum_{m=0}^{rN}\frac{\Phi_{\mathrm{bath}}^m}{m!}
=e^{\Phi_{\mathrm{bath}}}.
}
\tag{46}
$$

这是全部阶数的重求和，不是将杂化耦合截断到二阶。外代数总共只有 $2rN$ 个生成元，$\Phi_{\mathrm{bath}}$ 的次数为二，所以指数级数也在 $m=rN$ 后精确终止。

不同连续时间可以取到同一离散源，并不破坏证明。以上是源多项式的恒等式；某项因重复生成元而为零时，完整配对公式与指数展开中相应的多项式关系仍成立。整个推导从未把同一单元的源误当成彼此独立的连续时间变量。

### 5.5 提出分段常变量，得到 QuAPI 单元积分系数

现在才把式 (40) 中的分段常表示展开。由于 $a_\alpha(\tau)=a_{i\alpha}$ 在 $C_i$ 内固定，普通时间积分只作用于核，故

$$
\begin{aligned}
\Phi_{\mathrm{bath}}
&=\sum_{i,j=0}^{N-1}\sum_{\alpha\gamma}
\int_{C_i}d\tau\int_{C_j}d\tau'\,
\bar a_{i\alpha}\Delta_{\alpha\gamma}(\tau,\tau')a_{j\gamma}\\
&=\sum_{i,j=0}^{N-1}\sum_{\alpha\gamma}
\bar a_{i\alpha}\left[\int_{C_i}d\tau\int_{C_j}d\tau'\,
\Delta_{\alpha\gamma}(\tau,\tau')\right]a_{j\gamma}.
\end{aligned}
$$

这一步再次把有序指数中一直保留的离散源显式写出来。定义

$$
\boxed{
\Delta_{i\alpha,j\gamma}
=\int_{C_i}d\tau\int_{C_j}d\tau'\,
\Delta_{\alpha\gamma}(\tau,\tau'),
\qquad i,j=0,\ldots,N-1.
}
\tag{47}
$$

于是式 (46) 精确成为

$$
\boxed{\mathrm{IF}[\bar a,a]
=\exp\left[\sum_{i,j=0}^{N-1}\sum_{\alpha\gamma}
\bar a_{i\alpha}\Delta_{i\alpha,j\gamma}a_{j\gamma}\right].}
\tag{48}
$$

独立 bath 使 $\Delta_{i\alpha,j\gamma}=0$（$\alpha\ne\gamma$），因而 $\mathrm{IF}=\prod_{\alpha=1}^r\exp[\sum_{ij}\bar a_{i\alpha}\Delta_{i\alpha,j\alpha}a_{j\alpha}]$。各模可独立构造 IF；杂质模间相互作用仍完整保留在 $K_O$ 中。若去掉独立 bath 假设，保留式 (48) 的交叉模系数即可，不能继续按模相乘。

每个 $C_j$ 均为整宽 $\delta\tau$，总长严格为 $\beta$。式 (47) 就是 QuAPI 的核单元积分结构。逻辑顺序是：有限指数乘积的精确 Dyson 展开 → Gaussian 重求和 → 单元积分，而不是先定义一个独立连续源模型再断言它与乘积相同。

对一个 bath 能级 $\epsilon$，把其权重 $V_{\alpha k}V_{\gamma k}^*$ 暂时提出。用 $T_{\delta\tau}$ 表示三角区域积分，用 $A_{\delta\tau}$ 表示去掉单元间位移后的矩形积分；二者都有连续的零能量极限：

$$
T_{\delta\tau}(\epsilon)=\frac{e^{-\delta\tau\epsilon}-1+\delta\tau\epsilon}{\epsilon^2},
\quad T_{\delta\tau}(0)=\frac{(\delta\tau)^2}{2},
\qquad
A_{\delta\tau}(\epsilon)=\frac{2[\cosh(\delta\tau\epsilon)-1]}{\epsilon^2},
\quad A_{\delta\tau}(0)=(\delta\tau)^2.
\tag{49}
$$

按单元位置，式 (47) 给出

$$
\Delta^{(k)}_{i\alpha,j\gamma}=V_{\alpha k}V_{\gamma k}^*
\begin{cases}
(1-n_k)e^{-(i-j)\delta\tau\epsilon_k}A_{\delta\tau}(\epsilon_k),&i>j,\\
-n_ke^{(j-i)\delta\tau\epsilon_k}A_{\delta\tau}(\epsilon_k),&i<j,\\
(1-n_k)T_{\delta\tau}(\epsilon_k)-n_kT_{\delta\tau}(-\epsilon_k),&i=j.
\end{cases}
\tag{50}
$$

完整 $\Delta$ 对 $k$ 求和。跨单元时，时间顺序在整个矩形中固定，分别得到下三角正号和上三角负号；矩形积分的两个一维指数积分相乘即得 $A_{\delta\tau}$。同单元的晚时湮灭三角形给出 $(1-n)T_{\delta\tau}(\epsilon)$，早时湮灭三角形给出 $-nT_{\delta\tau}(-\epsilon)$，因此对角项是两项之和。

从式 (38) 还可直接检查热边界：将时间差从负值加上 $\beta$ 时，$(1-n_k)e^{-\beta\epsilon_k}=n_k$，故 $\Delta$ 的热延拓为反周期。这个 bath 收缩的热关系与第三节杂质迹的反周期闭合相容，但两者各有独立来源。

### 5.6 得到完整共同积分与真实物理目标的关系

把式 (33)、(48) 代入已经证明的式 (19)，得到

$$
\boxed{
\begin{aligned}
Z_{O,\mathrm{GTEMPO}}&=Z_{\mathrm{bath}}\int D[\bar a,a,\bar\chi,\chi]\,K_O
\exp\left[\sum_{i\alpha,j\gamma}
\bar a_{i\alpha}\Delta_{i\alpha,j\gamma}a_{j\gamma}\right],\\
Z_{\mathrm{GTEMPO}}&=Z_{\mathrm{bath}}\int D[\bar a,a,\bar\chi,\chi]\,K_I
\exp\left[\sum_{i\alpha,j\gamma}
\bar a_{i\alpha}\Delta_{i\alpha,j\gamma}a_{j\gamma}\right],\\
\langle O\rangle_\beta
&=\frac{\int D[\bar a,a,\bar\chi,\chi]\,K_O\mathrm{IF}[\bar a,a]}
{\int D[\bar a,a,\bar\chi,\chi]\,K_I\mathrm{IF}[\bar a,a]}+O(\delta\tau).
\end{aligned}
}
\tag{51}
$$

前两行是对 $U_{\mathrm{GTEMPO}}(\beta)$ 的精确表示，第三行才与真实热态比较。$Z_{\mathrm{bath}}$ 在比值中严格约去；不带该因子的积分对应约化配分函数 $Z_{\mathrm{GTEMPO}}/Z_{\mathrm{bath}}$，不是完整联合配分函数 $Z_{\mathrm{GTEMPO}}$。

以上已对任意固定有限自由正常态 bath、任意有限杂质模数、任意偶的含相互作用 $H_S$ 证明主线。若 bath 不守粒子数而有配对，Wick 仍可用，但须加入异常收缩；那是 Nambu 型 IF，不是本文式 (48) 的正常态 bath 接口。

## 6. 最后核对虚时间 QuAPI 与 GTEMPO 的实际目标

**本节的结论是：虚时间代码使用 $N$ 个整宽 IF 单元，没有额外的耦合传播端点；其路径闭合顺序与本文通过循环移动起点相联系，但当前 QuAPI 的小能量分支不完全等于解析单元积分。** 为逐项核验，先对照变量节点和边界，再分别检查 SIAM 精确系统核、通用短时核、QuAPI 系数和观测量插入。以下代码对应关系均指系数准确、张量表示不截断时的代数目标；谱积分求积与张量压缩误差另行计算。

首先核对两个关键对应关系：

$$
\boxed{
\begin{aligned}
a_j&\longleftrightarrow x_{j+2}\quad(j=0,\ldots,N-1),\\
U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta)
&=(\widetilde U_BU_S)^N
=U_S^{-1}U_{\mathrm{GTEMPO}}(\beta)U_S.
\end{aligned}
}
$$

$x_0,x_1$ 是只负责迹闭合与恒等分解的辅助节点，不进入 IF；上标 $\mathrm{cyc}$ 标明代码的乘积起点。下面先从源码确认这些关系，再核查系数与观测量插入是否也相同。

### 6.1 IF 的单元数及布局：虚时间没有额外的耦合传播端点

[虚时间格点定义](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/lattices/imaginarytime.jl:6) 令 `Nτ=N`、`k=kτ=N+1`、`β=N*δτ`，每个 band 存放编号 $0,1,\ldots,N+1$ 的 $N+2$ 对变量。但格点对数不等于耦合时间单元数。

[虚时间 IF 构造](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/influencefunctional/hybridization/partialif/imaginarytime.jl:6) 明确使用 `k=lattice.k-1=N`，遍历 `i=1:N`；调用 `partialif_hybrid(lattice,i+1,[0;view(corr,i,1:k)])`。因此实际 IF 为

$$
\mathrm{IF}_{\mathrm{pkg}}
=\exp\left[\sum_{i,j=0}^{N-1}
\bar x_{i+2}\Delta_{ij}x_{j+2}\right]
\tag{52}
$$

（先写单 band）。两个辅助节点 $x_0,x_1$ 均不进入 IF。对照表为：

| 角色              | 本文变量                                | package 节点             |
| --------------- | ----------------------------------- | ---------------------- |
| 迹辅助变量           | $\chi,\bar\chi$                     | $x_0,\bar x_0$         |
| 额外的杂质恒等分解变量     | 下一小节的 $\eta,\bar\eta$               | $x_1,\bar x_1$         |
| 第 $j$ 个整宽耦合时间单元 | $a_j,\bar a_j$                      | $x_{j+2},\bar x_{j+2}$ |
| 单元区间            | $C_j=[j\delta\tau,(j+1)\delta\tau)$ | QuAPI 的第 $j+1$ 行、列     |

这与实时文档中提到的额外的耦合传播端点不同：虚时间实现是 $N$ 个系统传播间隔、$N$ 个整宽耦合时间单元，总长均为 $\beta$。两个辅助变量不增加任何 $U_B$ 或 $\widetilde U_B$。源码构造中的 `vacuumstate(lattice)` 是 Grassmann 张量网络的常数多项式起点，不能据此把这里的配分函数理解成从物理杂质真空出发；热力学边界由式 (53) 的迹闭合给出。

### 6.2 系统核与迹边界：代码闭合的是哪个算符乘积

[边界构造](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/gvconnections/boundarycondition.jl:24) 添加 $e^{-\bar x_0x_{N+1}}e^{\bar x_1x_0}$；[系统连接](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/gvconnections/bulkconnection.jl:24) 的方向为从节点 $p$ 到 $p+1$。因此，用精确系统核时，代码的闭合系统因子为

$$
K_{\mathrm{pkg}}
=e^{-\bar x_0x_{N+1}+\bar x_1x_0}
\prod_{p=1}^{N}\langle\bar x_{p+1}|U_S|x_p\rangle.
\tag{53}
$$

多模时，边界指数对模求和，各系统核是完整多模核。将闭合路径的起点循环移动到系统传播之前，得到另一种乘积顺序；用上标 $\mathrm{cyc}$ 标明这种循环移动：

$$
U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta)=(\widetilde U_B U_S)^N.
\tag{54}
$$

为了说明代码多出的 $x_1$ 从何而来，从右往左代入每个 $\widetilde U_B$，并在最右端插入真正的恒等分解 $I=\int D[\bar\eta,\eta]|\eta\rangle\langle\bar\eta|$，得到

$$
\begin{aligned}
U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta)=\int D[\bar\eta,\eta]\prod_jD[\bar a_j,a_j]\,
|a_{N-1}\rangle&\left[\prod_{j=1}^{N-1}
\langle\bar a_j|U_S|a_{j-1}\rangle\right]\\
&\times\langle\bar a_0|U_S|\eta\rangle\langle\bar\eta|
U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}.
\end{aligned}
\tag{55}
$$

取迹，再用 $\chi$ 关闭路径，便得到式 (53)。这次 $x_1=\eta$ 是普通恒等分解的变量，$x_0=\chi$ 是迹变量；与所有 IF 的 Grassmann 变量保持独立。

因此，精确系数和不截断条件下，

$$
\boxed{
\int\prod_{p=0}^{N+1}D[\bar x_p,x_p]\,
K_{\mathrm{pkg}}\mathrm{IF}_{\mathrm{pkg}}
=\frac{\operatorname{Tr}U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta)}{Z_{\mathrm{bath}}}
=\frac{\operatorname{Tr}U_{\mathrm{GTEMPO}}(\beta)}{Z_{\mathrm{bath}}}.
}
\tag{56}
$$

最后一个等号在完整物理算符层面成立：$U_S$ 可逆，且 $U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta)=U_S^{-1}U_{\mathrm{GTEMPO}}(\beta)U_S$。这是相似变换后的普通矩阵迹恒等式，不是对含 Grassmann 奇系数的 bath 迹滥用循环性。

对一般观测量，有限步长下不能把这两个顺序直接混同：

$$
\operatorname{Tr}(O U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta))
=\operatorname{Tr}[(U_S O U_S^{-1})U_{\mathrm{GTEMPO}}(\beta)].
\tag{57}
$$

若要在代码布局中精确表示本文的 $Z_{O,\mathrm{GTEMPO}}=\operatorname{Tr}(O U_{\mathrm{GTEMPO}}(\beta))$，应将式 (53) 的边界重叠 $e^{\bar x_1x_0}$ 换为

$$
\boxed{
\langle\bar x_1|\widetilde O|x_0\rangle,
\qquad \widetilde O=U_S^{-1}OU_S.
}
\tag{58}
$$

若 $[O,H_S]=0$，则 $\widetilde O=O$，两种离散顺序在这个观测量上完全一致。否则，直接用 $O$ 边界表示的是 $\operatorname{Tr}(O U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta))$；由于 $U_S O U_S^{-1}=O+O(\delta\tau)$，它仍有 $O(\delta\tau)$ 的真实热态误差上界，但不是式 (51) 的同一个有限步长分子。

### 6.3 预定义 SIAM：这里的系统核是精确的

主定理使用精确 $U_S=e^{-\delta\tau H_S}$。它直接覆盖预定义 `AndersonIM`，其中源码参数对应

$$
H_S=\mu(n_\uparrow+n_\downarrow)+U_{\mathrm{int}}n_\uparrow n_\downarrow.
\tag{59}
$$

这里 $U_{\mathrm{int}}$ 对应源码字段 `model.U`，以免与传播子 $U_S$ 混淆；$\mu$ 是此处杂质能级参数，本文 bath 能量已按所选平衡态 Hamiltonian 定义，不能不作说明就再减一次化学势。

令 $q=e^{-\delta\tau\mu}$ 为单占据态的传播权重。直接在四个杂质 Fock 态上计算，可得精确核

$$
\begin{aligned}
\langle\bar y|e^{-\delta\tau H_S}|z\rangle
=\exp\big[&q(\bar y_\uparrow z_\uparrow+\bar y_\downarrow z_\downarrow)\\
&+q^2(e^{-\delta\tau U_{\mathrm{int}}}-1)
\bar y_\uparrow\bar y_\downarrow z_\downarrow z_\uparrow\big].
\end{aligned}
\tag{60}
$$

展开指数后，真空、单占据和双占据核的权重依次为 $1,q,q,q^2e^{-\delta\tau U_{\mathrm{int}}}$；其中两个二次项相乘本来已贡献 $q^2$ 的四次项，所以额外四次系数必须是 $q^2(e^{-\delta\tau U_{\mathrm{int}}}-1)$。

这正是 [SIAM 虚时间核](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/predefined/siam.jl:153) 及 [siam_coeffs](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/predefined/siam.jl:290) 的系数和四次项顺序。此接口无需再对 $U_S$ 做短时近似。多个独立正常 bath band 时各 band 的 IF 相乘；若 $\Delta$ 在 band 间非对角，则需要含交叉 band 系数的 IF，不能用独立单 band 构造冒充一般矩阵 bath。

### 6.4 通用 ImpurityHamiltonian：另有一个局部系统近似

[通用接口](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/generalimpurity/sysdynamics.jl:1) 使用参数 `dt=-δτ`，逐项乘上 Grassmann 项的指数，再由 [sysdynamics_util!](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/models/generalimpurity/sysdynamics.jl:70) 添加重叠连接。对其支持的正常序偶 Hamiltonian，未截断的单步核为

$$
\widehat G_S(\bar y,z)
=e^{\sum_\alpha\bar y_\alpha z_\alpha}
 e^{-\delta\tau H_S^{\mathrm{normal}}(\bar y,z)}.
\tag{61}
$$

这里 $H_S^{\mathrm{normal}}$ 是把正常序算符中的 $d^\dagger,d$ 按原顺序换成 $\bar y,z$ 的符号。它与本证明的 $Q$ 不是同一种映射。各 Grassmann 偶项彼此交换，因此逐项乘指数等于和的指数；这不是说物理 Hamiltonian 各项彼此交换。

设 $\widehat U_S$ 是核为 $\widehat G_S$ 的唯一杂质算符。其存在性与唯一性由 Fock 矩阵单位和核单项式的双射给出。正常序相干态核恒等式为

$$
\langle\bar y|H_S|z\rangle
=e^{\bar y z}H_S^{\mathrm{normal}}(\bar y,z).
\tag{62}
$$

逐个正常序单项式可用 $d|z\rangle=z|z\rangle$、$\langle\bar y|d^\dagger=\langle\bar y|\bar y$ 验证。展开式 (61) 的一次项，由有限维核映射的有界逆，得到

$$
\widehat U_S=I-\delta\tau H_S+O((\delta\tau)^2)=U_S(\delta\tau)+O((\delta\tau)^2).
\tag{63}
$$

例如单模 $H_S=\epsilon n$，通用核的系数是 $1-\delta\tau\epsilon$，精确核是 $e^{-\delta\tau\epsilon}$，两者一般不同。将 $U_S$ 换为 $\widehat U_S$ 后，式 (55) 的积分仍精确表示 $(\widetilde U_B\widehat U_S)^N$；式 (C2)–(C5) 再作一次局部替换估计即有全局 $O(\delta\tau)$ 上界。通用接口因此也能成立一阶理论，但不能把其系统核称作式 (60) 那样的精确系统核。

### 6.5 QuAPI 系数：解析区域吻合，但当前小能量分支有差异

[GTEMPO 的虚时间包装](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/correlationfunction.jl:8) 调用 `Δτ(bath,N=lattice.N,δτ=lattice.δτ)`。本证明要求 `bath.β=lattice.β=N*δτ`；代码在不相等时仅给警告，故调用者必须满足这个条件才能解释为同一个平衡态问题。

[QuAPI 的有效虚时间实现](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:53) 生成长度 $N$ 的两个数组。对默认色散、零 bath 化学势参数，并在未进入小能量分支时，将谱权重暂时提出，各项与式 (50) 对应为：

| 积分区域      | 数学结果                                                        | 当前有效函数位置                                                                                       |
| --------- | ----------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| $i>j$ 的矩形 | $(1-n)e^{-(i-j)\delta\tau\epsilon}A_{\delta\tau}(\epsilon)$ | [_fermionic_gⱼₖ_i](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:179) |
| $i<j$ 的矩形 | $-ne^{(j-i)\delta\tau\epsilon}A_{\delta\tau}(\epsilon)$     | [_fermionic_gₖⱼ_i](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:232) |
| 同单元晚时三角形  | $(1-n)T_{\delta\tau}(\epsilon)$                             | [_fermionic_gⱼⱼ_i](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:206) |
| 同单元早时三角形  | $-nT_{\delta\tau}(-\epsilon)$                               | [_fermionic_gₖₖ_i](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:261) |

`CorrelationMatrix` 的 [二维取值规则](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/corrmat.jl:18) 在非对角位置取相应数组元素，在对角位置返回 `ηⱼₖ[1]+ηₖⱼ[1]`，恰好是两个三角形之和。不能只取其中一个数组的首项。

但是，当前源码的小能量分支不能与式 (47) 全部严格认同。[QuAPI_tol](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/util.jl:12) 为 $10^{-6}$；在 `abs(ε)≤QuAPI_tol` 时，早时三角形函数使用了

```julia
# 当前有效代码，两条 ε 与 μ 的分支分别为：
-0.5*δτ^2 * x / 2 * (x+1)
-0.5*δτ^2 / 2 * (x+1)
```

位置分别是 [第 268 行](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:268) 和 [第 275 行](/home/zjsun/Desktop/Agent/try_sth/pkgs/QuAPI/src/imagtime/fermionic.jl:275)。按 Julia 通常的同优先级乘除从左往右规则，末尾乘 $(x+1)$ 不会成为分母。

以传入该核函数的 $\epsilon=\mu=0$ 为例，$x=1$、$n=1/2$，于是：

$$
\begin{array}{c|ccc}
&\text{晚时三角形}&\text{早时三角形}&\text{合并对角项}\\ \hline
\text{式 (47) 的精确值}&(\delta\tau)^2/4&-(\delta\tau)^2/4&0\\
\text{当前分支字面值}&(\delta\tau)^2/4&-(\delta\tau)^2/2&-(\delta\tau)^2/4
\end{array}
\tag{64}
$$

这已经在单能级核层面证明了差异。相应正确零能量极限应含 $x/(x+1)$ 或 $1/(x+1)$，而不是以上乘法写法。其他小能量分支也在用解析极限近似有限非零能量，不能把它们称作精确单元积分。

因此本证明的精确代码对应关系有明确范围：**节点、系统核、边界、矩形与三角形布局及非小能量解析系数可以闭合；当前小能量分支的字面输出与精确 QuAPI 单元积分不完全相等。** 对给定连续谱，这个分支对最终谱积分的影响还取决于谱权重和数值求积，本文不凭该核层面的反例推断任意实际计算的误差大小；也不把它吸收到已经证明的精确等号中。package 未作修改。

有限 bath 证明中只对能级权重求和。源码使用有界谱函数的数值积分接口；将有限和写成谱积分是系数表示的比较，不构成无限 bath 极限下的一致误差证明。默认之外的色散或化学势参数还必须与 Hamiltonian 中实际传播和占据的能量约定逐项一致，本文不作未经核验的全接口认同。

### 6.6 真实观测量与 Grassmann 变量的插入不是自动相同

式 (18) 的 $K_O$ 或式 (58) 是由真实算符 $O$ 推出的插入规则。[occupation 的虚时间接口](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/occupation.jl:15) 则通过 [Gτ 与 gf](/home/zjsun/Desktop/Agent/try_sth/pkgs/GTEMPO/src/observables/gf.jl:15) 计算归一化插入 $1+\bar x_i x_i$。必须区分节点角色。

单模在只作恒等分解的节点 $i=1$，$Q[1+\bar x_1x_1]=n$，故插入精确对应 $\operatorname{Tr}(nU_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta))$。迹辅助节点 $i=0$ 同样成立：对杂质矩阵单位验证

$$
\int D[\bar\chi,\chi](1+\bar\chi\chi)
\langle-\bar\chi|X|\chi\rangle=\operatorname{Tr}_S(nX).
\tag{65}
$$

若选择 参与 IF 的节点 $i=j+2$，插入实际将该单元的 $\widetilde U_B$ 换成

$$
\widetilde U_B^{[n]}(\delta\tau)=Q[(1+\bar a a)e^{-\delta\tau H_B^{\mathrm G}}],
\tag{66}
$$

它在有限 $\delta\tau$ 下一般不等于 $n\widetilde U_B$ 或 $\widetilde U_Bn$。但 $\widetilde U_B^{[n]}=n+O(\delta\tau)$、$n\widetilde U_B=n+O(\delta\tau)$，所以稳定性保证单次插入的分子差别为 $O(\delta\tau)$。

将插入处换成 $n\widetilde U_B U_S$ 后，其余因子都是 $\widetilde U_B U_S$ 的幂；在联合空间的迹中循环移动这些幂，便得到 $\operatorname{Tr}[n U_{\mathrm{GTEMPO}}^{\mathrm{cyc}}(\beta)]$。再用近似传播子的全局误差和同一分母估计，就得到真实占据数的 $O(\delta\tau)$ 上界。因此参与 IF 的节点的这个插入具有一阶上界，但不等于某个指定有限步长 $K_O$ 的精确算符插入。

对 SIAM 的 $n_\alpha$ 有 $[n_\alpha,H_S]=0$，故式 (57) 的乘积顺序问题不再造成额外差别；参与 IF 的节点上，Grassmann 插入与算符插入的差别仍需保留。多模的占据插入可对指定模逐模使用同一论证。任意偶观测量则应从其相干态算符核推导，不应仅凭变量标签猜测插入。

## 7. 推导结论与验证范围

**本节将前面的结果合成一条带误差标记的推导链，并界定验证支持的结论。** 先区分真实传播的近似与近似传播子的精确积分表示，再列出有限模型的代数核验。具有非零一阶误差的解析例子放在附录 F，用来区分“一阶上界”与“一般不是二阶”。

最终物理结论与代码核对条件分别为

$$
\boxed{
\langle O\rangle_\beta
=\frac{\int D[\bar a,a,\bar\chi,\chi]\,K_O\,\mathrm{IF}}
{\int D[\bar a,a,\bar\chi,\chi]\,K_I\,\mathrm{IF}}
+O(\delta\tau).
}
$$

其中积分精确表示本文的近似算符；要与 package 的有限步长结果认同，还须满足第六节的系统核、IF 系数及插入规则。下面将这些等号与近似号合并，并说明代数核验的范围。

### 7.1 已经闭合的逻辑链

对固定有限正常态自由 bath 和有限杂质，精确系统核的主结论可以写为

$$
\boxed{
\begin{aligned}
e^{-\beta H}
&=[U_S(\delta\tau)U_B(\delta\tau)]^N+O(\delta\tau)\\
&=[U_S(\delta\tau)\widetilde U_B(\delta\tau)]^N+O(\delta\tau),\\
\operatorname{Tr}(O U_{\mathrm{GTEMPO}}(\beta))
&=\int D[\bar a,a,\bar\chi,\chi]\,K_O\operatorname{Tr}_B(U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G})\\
&=Z_{\mathrm{bath}}\int D[\bar a,a,\bar\chi,\chi]\,K_O
 e^{\sum\bar a\Delta a}.
\end{aligned}
}
\tag{67}
$$

前两行包含两项实际近似，并已在使用积分表示之前证明其稳定性和全局上界；后两行全为恒等式。真实观测量的归一化误差由式 (14)、(C5) 控制。与当前代码对照时，还须应用第六节关于闭合路径起点和接口的条件，并单独保留式 (64) 的实际系数差异。

这里没有证明任意 bath 极限的一致收敛，没有证明任意张量截断策略随 $\delta\tau$ 自动收敛，也没有证明全部观测量都有非零的一阶误差系数。

### 7.2 代数核验及其限度

配套的 [verify_imaginary_order1.py](/home/zjsun/Desktop/Agent/try_sth/GTEMPO_theory/verify_imaginary_order1.py) 用 Python 有限矩阵与外代数作独立代数核验：检查虚时间积分近似的二次符号；比较指数乘积 bath 迹与二次指数的全部 Grassmann 系数；比较共同积分与近似传播子的取迹结果；核对一般非二阶的解析系数和式 (64) 的字面算术差异。

本次核验已通过：单模与双模积分近似的二次恒等式；13 组完整 IF 外代数比较（最高到八次，最大绝对差约 $2.22\times10^{-16}$）；9 组单模共同积分及代码布局比较；2 组含相互作用 SIAM 的共同积分比较；以及式 (F5) 的主误差系数和式 (64) 的字面算术反例。

这些是有限实例中的符号和实现审计，不能替代第二至第五节的一般证明，也不是实际运行 Julia package。源码位置以本次读取的文件版本为准。

## 附录 A. 相干态、积分测度与 Grassmann 多项式对应的算符

本附录固定正文积分的全部约定，并证明积分映射 $Q$ 的逐项对应表。先计算单模核，再逐模推广，最后说明 bath 奇算符系数和伴随的符号。

积分映射定义为 $Q[f]=\int D[\bar a,a]|a\rangle f\langle\bar a|$。单模取 $|1\rangle=d^\dagger|0\rangle$、$n=d^\dagger d$，并规定

$$
|a\rangle=e^{-a d^\dagger}|0\rangle,
\qquad\langle\bar a|=\langle0|e^{-d\bar a},
\qquad D[\bar a,a]=d\bar a\,da\,e^{-\bar a a},
\qquad\int d\bar a\,da\,a\bar a=1.
$$

使用分次张量积：交换不同因子的齐次对象产生 $(-1)^{pq}$，$p,q$ 为各自宇称。Grassmann 变量穿过奇 Fock 基矢产生负号。因此 $\langle\bar x|y\rangle=e^{\bar x y}$。Grassmann 单项式规范次序取 $a_1,\bar a_1,a_2,\bar a_2,\ldots$，bath 系数置最右；物理 Fock 次序为杂质模在前、bath 模在后。

不能把相干 ket、bra 的分次结构忽略后再套普通复数向量规则。一个安全的逐项检验是先计算核：

$$
\langle\bar x|Q[f]|y\rangle
=\int d\bar a\,da\,(1+a\bar a)(1+\bar x a)
f(\bar a,a)(1+\bar a y).
$$

下表的每项都由将 $f$ 代入后抽取 $a\bar a$ 的系数得到：

| $f$        | 积分所得核        | 对应算符        |
| ---------- | ------------ | ----------- |
| $1$        | $1+\bar x y$ | $I$         |
| $a$        | $y$          | $d$         |
| $\bar a$   | $\bar x$     | $d^\dagger$ |
| $a\bar a$  | $1$          | $I-n$       |
| $\bar a a$ | $-1$         | $n-I$       |

例如 $f=1$ 的两个贡献为测度给出的 $1$，以及两个重叠因子的交叉项给出的 $\bar x y$；$f=a\bar a$ 已占满积分生成元，其他非恒定项全为零。对于 $f=a$，唯一能占满积分变量的项是 $a\bar a y$，得到 $y$；对于 $f=\bar a$，对应项是 $\bar x a\bar a$，得到 $\bar x$。最后一行由 $\bar a a=-a\bar a$ 得到。任一单模算符的四个 Fock 矩阵元分别对应核中四个独立单项式，故核相等唯一确定算符相等，不只是比较少数期望值。

多模定义为上述相干态的有序分次张量积，测度为

$$
D[\bar a,a]=\prod_{\alpha=1}^{r}d\bar a_\alpha\,da_\alpha\,
e^{-\sum_\alpha\bar a_\alpha a_\alpha}.
$$

对按模次序排列的单项式逐模积分，每一步就是上表，得到同次序的 $I,d_\alpha,d_\alpha^\dagger,I-n_\alpha$ 的分次乘积。测度对为偶；将一个模块移过其余模时，积分端和算符端的奇因子交换符号相同。这给出按模数归纳的推广，并给出

$$
Q[a_\alpha\bar a_\gamma]=d_\alpha d_\gamma^\dagger,
\qquad
Q[\bar a_\alpha a_\gamma]
=d_\alpha^\dagger d_\gamma-\delta_{\alpha\gamma}I.
$$

多模核的独立单项式同样逐一标记全部 Fock 矩阵元，故核映射仍为双射。含 bath 奇系数时，先移到最右，再用此表：例如 $F^\dagger a=-aF^\dagger$，所以 $Q[F^\dagger a]=-dF^\dagger=F^\dagger d$。Grassmann 单项式 $m$ 与最右系数 $C$ 的推广为 $Q[mC]=Q[m]C$，这里的乘积是物理 CAR 代数中的乘积。若用普通矩阵张量积，bath 奇算符含完整杂质宇称串，绝不能省略。

取 Grassmann 变量的伴随时反转乘积次序，并取 $a^\dagger=\bar a$。上表逐项保持伴随；多模反转后重新排序在两端产生相同符号，最右 bath 系数也按分次规则处理，故 $Q[f^\dagger]=Q[f]^\dagger$。本证明不把虚时间的相互作用绘景与伴随互换。

$Q$ 线性但不保乘法。例如 $Q[\bar a a]=n-I\ne Q[\bar a]Q[a]=n$。特别是 $\widetilde U_B(\delta\tau)$ 一般不等于 $e^{-\delta\tau Q[H_B^{\mathrm G}]}$。

## 附录 B. 一阶分裂的局部误差、乘积稳定性与归一化

本附录证明式 (7)、(8)。所有范数均为有限联合 Fock 空间上的算符范数，令 $M=\|H_S\|+\|H_B\|$。由矩阵指数级数和范数的次乘性，对 $s\geq0$ 有

$$
\|e^{-sH}\|\leq e^{sM},\qquad
\|e^{-sH_S}e^{-sH_B}\|\leq e^{sM},\qquad
\|(e^{-\delta\tau H})^k\|,\ \|(U_S U_B)^k\|\leq e^{k\delta\tau M}.
\tag{B1}
$$

因而当 $k\delta\tau\leq\beta$ 时，任意这些部分乘积的范数均不超过 $e^{\beta M}$。这个上界与步数无关，是本节需要的虚时间乘积稳定性。

为直接证明单步误差，引入连接两个单步传播子的辅助函数 $\Phi_{\mathrm{split}}(s)=e^{-(\delta\tau-s)H}e^{-sH_S}e^{-sH_B}$，$0\leq s\leq \delta\tau$。其端点分别为 $\Phi_{\mathrm{split}}(0)=e^{-\delta\tau H}$、$\Phi_{\mathrm{split}}(\delta\tau)=U_S U_B$。对它求导，得到

$$
\frac{d\Phi_{\mathrm{split}}(s)}{ds}
=e^{-(\delta\tau-s)H}[H_B,e^{-sH_S}]e^{-sH_B}.
\tag{B2}
$$

其中的交换子也可通过对两个指数之间插入 $H_B$ 的表达式求导、再积分而精确写成

$$
[H_B,e^{-sH_S}]
=\int_0^s du\,
e^{-(s-u)H_S}[H_S,H_B]e^{-uH_S}.
\tag{B3}
$$

由式 (B2)、(B3) 及各指数的范数上界，$\|\Phi_{\mathrm{split}}'(s)\|\leq s e^{\delta\tau M}\|[H_S,H_B]\|$。在 $0\leq s\leq \delta\tau$ 上积分，得到单步估计

$$
\|e^{-\delta\tau H}-U_S U_B\|
\leq\frac{(\delta\tau)^2}{2}e^{\delta\tau M}\|[H_S,H_B]\|.
\tag{B4}
$$

这证明式 (5)，没有使用任何传播子的幺正性。接着用非交换乘积的伸缩恒等式

$$
(e^{-\delta\tau H})^N-(U_S U_B)^N
=\sum_{j=0}^{N-1}
(e^{-\delta\tau H})^{N-1-j}[e^{-\delta\tau H}-U_S U_B](U_S U_B)^j.
\tag{B5}
$$

每个求和项的范数至多为 $\frac{(\delta\tau)^2}{2}e^{N\delta\tau M}\|[H_S,H_B]\|$。由于 $N\delta\tau=\beta$，故

$$
\|e^{-\beta H}-U_{\mathrm{split}}(\beta)\|
\leq\frac{\beta \delta\tau}{2}e^{\beta M}\|[H_S,H_B]\|
=:\varepsilon_{\mathrm{split}}=O(\delta\tau).
\tag{B6}
$$

最后将这个界传递到两个迹及其比值。设联合 Fock 空间维数为 $D$。由 $|\operatorname{Tr}A|\leq D\|A\|$，有

$$
|\operatorname{Tr}U_{\mathrm{split}}(\beta)-Z|\leq D\varepsilon_{\mathrm{split}},
\qquad
|\operatorname{Tr}[O U_{\mathrm{split}}(\beta)]-Z_O|
\leq D\|O\|\varepsilon_{\mathrm{split}}.
\tag{B7}
$$

固定模型及 $\beta$ 时，真实的 $Z>0$ 与 $\delta\tau$ 无关。因此，对充分小的 $\delta\tau$，$|\operatorname{Tr}U_{\mathrm{split}}(\beta)|\geq Z-D\varepsilon_{\mathrm{split}}\geq Z/2$。利用真实热态的 $|\langle O\rangle_\beta|\leq\|O\|$，得到

$$
\left|\frac{\operatorname{Tr}[O U_{\mathrm{split}}(\beta)]}{\operatorname{Tr}U_{\mathrm{split}}(\beta)}-\frac{Z_O}{Z}\right|
\leq
\frac{|\operatorname{Tr}[O U_{\mathrm{split}}(\beta)]-Z_O|
+\|O\|\,|\operatorname{Tr}U_{\mathrm{split}}(\beta)-Z|}{|\operatorname{Tr}U_{\mathrm{split}}(\beta)|}
\leq\frac{4D\|O\|}{Z}\varepsilon_{\mathrm{split}}=O(\delta\tau).
\tag{B8}
$$

以上证明的是误差上界，不能据此断言某个迹或某个观测量的首个非零误差必定正比于 $\delta\tau$。这里的常数针对固定有限 bath；尚未对浴极限作一致性结论。

## 附录 C. 耦合传播近似的局部误差与全局误差

本附录证明式 (12)、(14)：先核验二次差额，再控制 $N$ 步乘积的范数，最后将误差传递到配分函数、观测量分子及其比值。

对 Grassmann 多项式各规范单项式的 bath 系数取算符范数之和，得到次乘性的有限维代数范数；$Q$ 为有界线性映射。含 Grassmann 变量的指数 Taylor 级数绝对收敛，余项可估为常数乘 $(\delta\tau)^3$，故式 (12) 是算符范数意义的估计。

为核验二次项，记含 Grassmann 变量的杂化为 $H_{\mathrm{hyb}}^{\mathrm G}=\sum_\alpha(\bar a_\alpha F_\alpha+F_\alpha^\dagger a_\alpha)$。含 $H_{\mathrm{bath}}$ 的项线性积分后完全相同。杂化平方中，两个产生算符对应的项、两个湮灭算符对应的项、以及“湮灭在左、产生在右”的项也相同；差别来自

$$
(\bar a_\alpha F_\alpha)(F_\gamma^\dagger a_\gamma)
=\bar a_\alpha a_\gamma F_\alpha F_\gamma^\dagger,
\qquad
(d_\alpha^\dagger F_\alpha)(F_\gamma^\dagger d_\gamma)
=d_\alpha^\dagger d_\gamma F_\alpha F_\gamma^\dagger.
$$

因此，由附录 A 中的 $Q[\bar a_\alpha a_\gamma]$，

$$
\boxed{
\begin{aligned}
Q[(H_B^{\mathrm G})^2]
&=H_B^2-\sum_\alpha F_\alpha F_\alpha^\dagger,\\
\widetilde U_B(\delta\tau)-U_B(\delta\tau)
&=-\frac{(\delta\tau)^2}{2}\sum_\alpha F_\alpha F_\alpha^\dagger
+O((\delta\tau)^3).
\end{aligned}
}
\tag{C1}
$$

单模也可直接检查：$(H_{\mathrm{hyb}}^{\mathrm G})^2=\bar a a FF^\dagger+a\bar a F^\dagger F$，对应于 $(n-I)FF^\dagger+(I-n)F^\dagger F$，而 $H_{\mathrm{hyb}}^2=nFF^\dagger+(I-n)F^\dagger F$。这个负号来自虚时间指数的二次系数为 $+(\delta\tau)^2/2$。

取固定 $\delta\tau_{\max}>0$ 足够小，使 $0<\delta\tau\leq \delta\tau_{\max}$ 时 $\|\widetilde U_B-U_B\|\leq c(\delta\tau)^2$。沿用附录 B 的 $M=\|H_S\|+\|H_B\|$，不再给单步乘积另取名字。直接有

$$
\|(U_S\widetilde U_B)-(U_SU_B)\|\leq ce^{\delta\tau\|H_S\|}(\delta\tau)^2,
\qquad
\|(U_S\widetilde U_B)\|\leq e^{\delta\tau M}(1+c(\delta\tau)^2)
\leq e^{\delta\tau M+c(\delta\tau)^2}.
\tag{C2}
$$

于是 $\|(U_S\widetilde U_B)^k\|\leq e^{k\delta\tau M+ck(\delta\tau)^2}$，$k\delta\tau\leq\beta$ 时一致有界；并未假定 $\widetilde U_B$ 是收缩或幺正算符。伸缩求和给出

$$
\|U_{\mathrm{GTEMPO}}(\beta)-U_{\mathrm{split}}(\beta)\|
\leq c\beta \delta\tau\,e^{\beta M+c\beta \delta\tau}.
\tag{C3}
$$

与式 (B6) 相加，记

$$
\varepsilon_{\mathrm{total}}=\frac{\beta \delta\tau}{2}e^{\beta M}\|[H_S,H_B]\|
+c\beta \delta\tau e^{\beta M+c\beta \delta\tau},
\qquad \|U_{\mathrm{GTEMPO}}(\beta)-e^{-\beta H}\|\leq\varepsilon_{\mathrm{total}}=O(\delta\tau).
\tag{C4}
$$

联合空间维数仍为 $D$，故 $|Z_{\mathrm{GTEMPO}}-Z|\leq D\varepsilon_{\mathrm{total}}$、$|Z_{O,\mathrm{GTEMPO}}-Z_O|\leq D\|O\|\varepsilon_{\mathrm{total}}$。当 $D\varepsilon_{\mathrm{total}}\leq Z/2$ 时，$|Z_{\mathrm{GTEMPO}}|\geq Z/2$，同式 (B8) 得

$$
\left|\frac{Z_{O,\mathrm{GTEMPO}}}{Z_{\mathrm{GTEMPO}}}-\frac{Z_O}{Z}\right|
\leq\frac{4D\|O\|}{Z}\varepsilon_{\mathrm{total}}.
\tag{C5}
$$

这些界不要求离散算符在任意粗步长下保持正性。它们针对固定有限 bath，不能直接当成 bath 模数趋于无穷时的一致估计。

## 附录 D. 因子重排与相干态求迹的符号

本附录证明第三节等式链使用的两条规则：含 Grassmann 变量的 bath 传播子可以穿过杂质偶对象，但不能交换彼此的次序；普通迹的相干态闭合必须保留负号。

### D.1 保持 bath 算符的乘积顺序

$U_{B,j}^{\mathrm G}$ 每一项的 Grassmann 系数宇称与 bath 系数宇称相同；相干 ket、bra 作为包含其 Grassmann 系数的整体也为偶。在分次张量积中，不同物理因子上的这些整体可以交换：交换 Grassmann 变量、交换物理对象及其交叉交换的总符号为 $+1$。此规则也允许 $U_{B,j}^{\mathrm G}$ 穿过偶杂质算符、系统核。若有重复生成元，相关单项式在两边都为零。

不同 $U_{B,j}^{\mathrm G}$ 的 bath 系数却作用在同一空间，通常不交换；不得改变 $U_{B,N-1}^{\mathrm G}\cdots U_{B,0}^{\mathrm G}$ 的次序。这里也没有对含 Grassmann 奇系数的 bath 迹使用普通循环性。

### D.2 为什么普通迹的闭合处有负号

正文中使用的相干态求迹公式是

$$
\operatorname{Tr}_S X
=\int D[\bar\chi,\chi]\,\langle-\bar\chi|X|\chi\rangle.
\tag{D1}
$$

单模时，分别代入矩阵单位 $|0\rangle\langle0|$、$|1\rangle\langle1|$，核为 $1$、$-\bar\chi\chi$；带测度积分均为 $1$。两个非对角矩阵单位的积分为零。多模逐个矩阵单位使用相同计算，得到对角元之和。由线性性式 (D1) 成立；带独立外部 Grassmann 变量时，非对角项仍积分为零，对角项的偶系数可按固定次序提出。

删除 bra 的负号会让占据态贡献变成 $-1$，所得是宇称加权迹，不是热力学配分函数。

## 附录 E. 有限自由 bath 的 Wick 递推证明

本附录为第五节提供全阶 Wick 定理的有限 bath 证明。先将一个 bath 湮灭算符移到乘积末端，再用自由热权重将末端项移回左侧，由此得到可归纳的收缩递推。

### E.1 纯 bath 线性算符的 Wick 递推

这一步不要求连续浴。设 $X_2,\ldots,X_{2m}$ 均为 $b_k,b_k^\dagger$ 的线性组合。将 $b_k$ 用 CAR 逐个移到最右，得到各反对易子的和，另有负的末端项。此时对纯 bath 矩阵才使用普通迹循环性。由 $b_k\rho_{\mathrm{bath}}^{\mathrm{free}}=e^{-\beta\epsilon_k}\rho_{\mathrm{bath}}^{\mathrm{free}} b_k$，末端期望等于 $e^{-\beta\epsilon_k}\langle b_k X_2\cdots X_{2m}\rangle_{\mathrm{bath}}$。移项得

$$
\langle b_k X_2\cdots X_{2m}\rangle_{\mathrm{bath}}
=(1-n_k)\sum_{j=2}^{2m}(-1)^j
\{b_k,X_j\}\,
\langle X_2\cdots\widehat X_j\cdots X_{2m}\rangle_{\mathrm{bath}}.
$$

帽号表示删掉该因子。以 $b_k^\dagger$ 开头的递推同理，只需将 $1-n_k$ 换成 $n_k$。因此对任一线性 $X_1$，递推系数恰好是 $\langle X_1X_j\rangle_{\mathrm{bath}}$。按 $m$ 归纳就得到全阶 Wick 配对公式。奇阶由宇称为零。这里迹循环时已经没有外部 Grassmann 变量，未假设含 Grassmann 变量的 bath 迹可普通循环。

### E.2 为什么完整杂化因子的配对公式全部写成加号

本小节补足第五节从具体四因子计算到一般配对公式的分次代数。先固定 $\tau_1>\cdots>\tau_{2m}$，将每个 $H_{\mathrm{hyb}}^{\mathrm G}(\tau_i)$ 按杂质模及 bath 算符类型逐项展开。每一项均写成 $\theta_iX_i$，源在左、bath 奇算符在右：对于 $\bar a_\alpha F_\alpha$，取 $\theta_i=\bar a_\alpha$、$X_i=F_\alpha$；对于 $F_\alpha^\dagger a_\alpha$，取 $\theta_i=-a_\alpha$、$X_i=F_\alpha^\dagger$。这些记号只在本附录用于证明符号，不改变正文的因子定义。

从乘积中提出全部源，共需 $(2m)(2m-1)/2$ 次奇对象交换，所以

$$
\left\langle(\theta_1X_1)\cdots(\theta_{2m}X_{2m})\right\rangle_{\mathrm{bath}}
=(-1)^m\theta_1\cdots\theta_{2m}
\left\langle X_1\cdots X_{2m}\right\rangle_{\mathrm{bath}}.
\tag{E1}
$$

对一个配对 $\mathcal P=\{(p_1,q_1),\ldots,(p_m,q_m)\}$，每对取 $p_k<q_k$；记 $\operatorname{sgn}(\mathcal P)$ 为将原标签次序排成 $p_1,q_1,\ldots,p_m,q_m$ 的置换符号。纯 bath Wick 项为 $\operatorname{sgn}(\mathcal P)\prod_k\langle X_{p_k}X_{q_k}\rangle_{\mathrm{bath}}$。另一方面，每对的完整平均为

$$
\left\langle(\theta_pX_p)(\theta_qX_q)\right\rangle_{\mathrm{bath}}
=-\theta_p\theta_q\left\langle X_pX_q\right\rangle_{\mathrm{bath}}.
$$

将各对相乘并把源恢复到原次序，恰有

$$
\begin{aligned}
&\prod_{k=1}^{m}
\left\langle(\theta_{p_k}X_{p_k})(\theta_{q_k}X_{q_k})\right\rangle_{\mathrm{bath}}\\
&\quad=(-1)^m\operatorname{sgn}(\mathcal P)
\theta_1\cdots\theta_{2m}
\prod_{k=1}^{m}\left\langle X_{p_k}X_{q_k}\right\rangle_{\mathrm{bath}}.
\end{aligned}
\tag{E2}
$$

它与式 (E1) 中该 Wick 配对的贡献完全一致。因此对所有配对求和时，使用完整二因子平均之积，每一项前面都写加号。再对所有算符类型、杂质模指标求和，即得正文的配对公式。二因子平均是 Grassmann 偶标量，配对块之间可交换；这并不意味着原 bath 算符可交换。

上述证明在每个严格时间排序区域成立。对任意时间标签，先用 $\mathcal T$ 排序整个偶因子；每对内部也按同一时间排序，因而得到全时间立方体上的同一公式。等时边界不影响积分。所有步骤都是外代数恒等式，即使不同 $\theta_i$ 是同一个生成元也成立。

## 附录 F. 一个具有非零一阶误差的虚时间例子

为了区分“一阶上界”与“一般并非二阶”，选取一个可在有限维中解析计算的例子：一个 impurity mode，一个零能量 bath mode，$H_S=H_{\mathrm{bath}}=0$，实数耦合 $v\ne0$，即 $F=vb$。于是

$$
H=H_B=v(d^\dagger b+b^\dagger d),\qquad
U_{\mathrm{GTEMPO}}(\beta)=\widetilde U_B(\delta\tau)^N.
\tag{F1}
$$

这里没有 Lie 分裂误差，也没有多余端点。由式 (C1) 对近单位矩阵取对数，

$$
\begin{aligned}
\log\widetilde U_B(\delta\tau)
&=-\delta\tau H_B-\frac{(\delta\tau)^2v^2}{2}(1-b^\dagger b)
+O((\delta\tau)^3),\\
\widetilde U_B(\delta\tau)^N
&=\exp\left[-\beta H_B-\frac{\beta\delta\tau v^2}{2}(1-b^\dagger b)
+O((\delta\tau)^2)\right].
\end{aligned}
\tag{F2}
$$

矩阵对数在 $\delta\tau$ 足够小时解析，余项乘 $N$ 后是 $O((\delta\tau)^2)$。对指数作一阶 Duhamel 展开并取普通算符迹，各插入位置经迹循环后相同，因此

$$
Z_{\mathrm{GTEMPO}}=Z-\frac{\beta\delta\tau v^2}{2}\operatorname{Tr}[e^{-\beta H_B}(1-b^\dagger b)]+O((\delta\tau)^2).
\tag{F3}
$$

在 Fock 次序 $|00\rangle,|10\rangle,|01\rangle,|11\rangle$ 下，$H_B$ 仅在两个单占据态之间交换，故

$$
Z=2+2\cosh(\beta v),\qquad
\operatorname{Tr}[e^{-\beta H_B}(1-b^\dagger b)]
=1+\cosh(\beta v)=Z/2.
\tag{F4}
$$

从而

$$
\boxed{\frac{Z_{\mathrm{GTEMPO}}-Z}{Z}=-\frac{\beta v^2}{4}\delta\tau+O((\delta\tau)^2).}
\tag{F5}
$$

这证明本文的 Grassmann 积分近似在一般情形下不能称为二阶；不是把实时半宽端点的误差系数搬到虚时间。该例子的精确 IF 对角系数为零，按式 (50) 的解析零能量极限定义；**它检验的是数学离散化，不是当前 QuAPI 错误小能量分支的输出。** 某些归一化观测量仍可能因对称性消掉一阶项，式 (F5) 不声称它们全部只有一阶。
