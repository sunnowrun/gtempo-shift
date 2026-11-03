# 基于 GTEMPO 理论的顶点移动复核

结论与证据见 [REPORT.md](REPORT.md)，完整数学链见 [derivation.md](docs/derivation.md)。本目录独立保存本次复核；原工作区报告和原始结果未改写。

## 环境与范围

数值计算在 CPU 服务器 `/public/home/zjsun/Agent-workspace/gc-ask/gtempo-shift` 进行。项目使用该目录已有的 `env/Project.toml`、`env/Manifest.toml` 和 `vendor/` 源码。作者六个仓库的提交见 [source_audit.json](results/source_audit.json)。GTEMPO 与 QuAPI HEAD 在本次复核时再次只读检查，未变化。注册依赖使用已锁定版本，未声称全是最新版本。

Python 使用 NumPy、SciPy；Julia 1.12.6 使用 GTEMPO、QuAPI、ImpurityModelBase、Z2Tensors 等锁定依赖。所有运行关闭本目录外的 Python 字节码写入。两个 `reference/` 工具从既有理论 fork 逐字复制，哈希核验通过，未修改外部文件。

本次区分三类验证：有限 Fock 矩阵的理论审计；完整非相互作用 Gaussian 轮廓积分；真实最新版 Julia GMPS 构造和收缩。`package.jl` 的实轴和 mixed 验证将显式有序复时间单元存放在 ImagGrassmannLattice 中，复用 GMPS 积分原语，**并非**声称调用了原生 mixed TTI 接口。`native.jl` 与 `alternatives.jl` 则直接使用真实实轴/虚轴/mixed 布局、闭合与观测接口。

## 复现

在服务器工作区根目录启动一个持久 Julia REPL，整个流程复用该进程。按站点规则申请 CPU 资源；以下环境路径均在授权服务器工作目录内：

```bash
export TMPDIR=/public/home/zjsun/Agent-workspace/gc-ask/gtempo-shift/tmp
export JULIA_DEPOT_PATH=/public/home/zjsun/Agent-workspace/gc-ask/gtempo-shift/depot:/public/home/zjsun/.julia
export JULIA_HISTORY=/public/home/zjsun/Agent-workspace/gc-ask/gtempo-shift/theory_reassessment/logs/julia_history.jl
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
/public/home/zjsun/julia-1.12.6/bin/julia --project=env --startup-file=no
```

Julia 中先导出最新版真实系数，再生成 Python 目标数据：

```julia
using GTEMPO, QuAPI, LinearAlgebra
BLAS.set_num_threads(1)
include("theory_reassessment/scripts/export.jl")
TheoryExport.run()
```

同一计算节点上执行 Python 数值脚本；先后依赖如下：

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 theory_reassessment/scripts/foundation.py
python3 theory_reassessment/scripts/operator_audit.py
python3 theory_reassessment/scripts/contour_audit.py
python3 theory_reassessment/scripts/package_inputs.py
python3 theory_reassessment/scripts/alternative_inputs.py
```

回到同一个 Julia 会话，继续：

```julia
include("scripts/study.jl")
include("scripts/kernel_audit.jl")
include("theory_reassessment/scripts/package.jl")
include("theory_reassessment/scripts/native.jl")
include("theory_reassessment/scripts/alternatives.jl")
include("theory_reassessment/scripts/benchmark.jl")
TheoryPackage.run();
TheoryNative.run();
TheoryAlternatives.run();
TheoryBenchmark.run();
```

`study.jl` 为原研究目录的原始脚本，只在进程内覆盖相关读取方法，以显式选择 raw 或旧 vertex 对照，不修改包源码。TheoryPackage 使用自身类型分派和显式矩阵，TheoryAlternatives 使用逐行 GMPS 原语，不依赖这个覆盖来构造新方案。TheoryNative 结束后恢复原研究的 `MODE=:current`。不要在生产 Julia 会话中混用这套研究覆盖。

最后在服务器执行：

```bash
python3 theory_reassessment/scripts/validate.py
```

本地仅作外部文档哈希和 Markdown 检查：

```bash
python3 theory_reassessment/scripts/source_audit.py
python3 theory_reassessment/scripts/check_math.py
```

重复执行数值脚本会覆盖本目录同名结果，应先保存已有结果。本次实际执行日志、失败控制和会话信息见 [execution.md](docs/execution.md)。大网格 Python 审计不包含 MPS 截断；实际 MPS 小链使用最大键维 160、截断阈值 10⁻¹³。

## 数据索引

| 文件 | 验证对象 |
|---|---|
| `cases.csv`、每例节点及 raw/vertex 矩阵 | 上游实际核和显式节点 |
| `source_blocks.csv` | 有限源指数与上游各块逐项比较 |
| `local_algebra.csv` | Q、正常链接、Weyl 局域恒等式及阶数 |
| `leading_error.csv` | 占据一阶系数、主体与端点分解 |
| `green_leading.csv` | 原大于 Green 的独立算符值与逐节点主误差 |
| `operator_convergence.csv` | 单模、相互作用双模的三轮廓全局收敛 |
| `link_and_sign.csv` | 一致链接、半宽端点、旧符号错误的判别实验 |
| `native.csv` | 上游原生网络的 raw/vertex 收缩 |
| `package.csv` | 34 个新网络与独立 Fock 目标的分子分母比较 |
| `alternatives.csv` | 半宽和一致链接在原生布局中的收缩 |
| `warm_timings.csv` | 同一 Julia 进程中预热后三次时间中位数 |
| `validation.json` | 自动断言及汇总，误差定义由脚本固定 |

表中的 NumPy 浮点数是确定性计算结果，没有统计误差条。收敛阶取最细两点比值；偶然小误差常数不能代替一般阶数证明。两种方案的不同观测量、不同时间长度、不同浴参数测试不得混成一条同成本排名。
