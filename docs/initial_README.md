# GTEMPO / QuAPI 顶点移动分析

主要结论和结果表见 [REPORT.md](REPORT.md)，推导见 [docs/mathematics.md](docs/mathematics.md)，图见 [results/convergence.png](results/convergence.png)。

## 工作目录与可复现记录

本地工作目录为 `/home/zjsun/Desktop/Agent/try_sth/test/gc-ask/gtempo-shift`；CPU 服务器目录为 `/public/home/zjsun/Agent-workspace/gc-ask/gtempo-shift`。全部数值计算在 CPU 服务器进行。对工作目录外只读取了 SSH 配置、Julia 安装和已有包缓存；新的包环境、编译缓存、临时文件、日志均写入指定工作目录。

测试时最新作者仓库源码在 `vendor/`，提交记录见 `results/source_versions.txt`，注册依赖版本见 `env/Manifest.toml`。原包源码未修改；研究性修改仅发生在持续 Julia 进程中，并在对照结束后恢复。

首次复现，在服务器工作目录内启动有计算资源的 Julia 会话，执行 `include("scripts/setup.jl")`。此脚本重新解析本工作目录内的路径依赖；已有 Manifest 中路径来自原服务器，不宜原样移植到其他机器。随后执行 `include("scripts/run_all.jl")`。这会重写同名输出，保留本次证据时请先复制到另一个结果目录。

绘图和汇总在同一 CPU 计算资源内运行 `python3 scripts/summarize.py`。Matplotlib 缓存目录应设在本工作目录内，例如设置 `MPLCONFIGDIR` 为本工作目录的 `mplconfig` 子目录。所有记录都具有确定参数；没有随机物理采样，随机系数抽查使用固定种子。

## 脚本

| 文件 | 目的 |
|---|---|
| `scripts/study.jl` | 三轮廓 PartialIF 收敛、独立精确参考、Richardson、真实浴尾部 Prony |
| `scripts/kernel_audit.jl` | 核耦合标度、从 MPS 逐项读取核、当前 TTI 对照 |
| `scripts/delay_mpo.jl` | 不含小极点除法的延迟状态核生成器 |
| `scripts/link_if.jl` | 精确单步 IF、匹配系统相位、对称链及失败控制 |
| `scripts/repaired_tti.jl` | 体尾部修正、精确局部和边界补偿、系数抽查 |
| `scripts/summarize.py` | 原始 CSV 的收敛阶与图，无平滑和选择性删点 |
| `scripts/check_math.py` | 检查本次 Markdown 编辑的数学定界符 |

## 会话交接

Slurm 作业 `313869`，节点 `node1`，4 CPU、16 GiB；分区默认和上限都是 6 天，本次未另设时限。Julia OS PID `113415`，终端工具会话 `62467`，项目为服务器的 `env`。会话在测试完成后保持存活，等待输入；模块已加载。改变已加载结构体或依赖版本时需要重启，普通参数扫描无需重启。

正式数值日志在服务器 `logs/`，本地同步副本为 `logs/`。日志包含首次编译、研究性修复和失败控制；失败探索没有当作通过的生产实现。核心结果足以由磁盘脚本重建，无需依赖会话里的唯一状态。
