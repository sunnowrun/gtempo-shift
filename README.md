# gtempo-shift

基于 GTEMPO 一阶理论，对 QuAPI / GTEMPO 顶点移动、混合轮廓、Prony 兼容性及 Weyl/K-contact 替代方案的数学与数值复核。

**最新结论与完整推导： [理论复核报告](theory_reassessment/REPORT.md) · [详细推导](theory_reassessment/docs/derivation.md) · [验证汇总](theory_reassessment/results/validation.json)**

实轴移动有明确的离散含义，能显著改善本次基准的一阶误差常数，但不一般升阶。旧 mixed 移动的附加符号令正向支杂化反号，趋向错误的连续轮廓。跨轮廓统一二阶路线采用原物理 IF、Weyl/K-contact 重建与精确系统半步；最小改动方案是 raw 结果的两网格外推。实际实现的范围、性能代价和未完成的 TTI 接口适配，均在报告中说明。

## 项目内容

| 路径 | 内容 |
|---|---|
| [theory_reassessment/](theory_reassessment/README.md) | 最终理论复核、34 个新 GMPS 网络、12 个原生布局对照及完整数据 |
| [REPORT.md](REPORT.md) 与 [docs/](docs/) | 第一阶段报告和探索记录；其中混合轴的解释由最终报告更新 |
| `scripts/`、`results/`、`logs/` | 第一阶段数值脚本、原始数据、图表与执行日志 |
| `env/`、[requirements.txt](requirements.txt) | Julia 环境锁定与 Python 版本 |
| [vendor/](THIRD_PARTY.md) | 六个作者仓库的完整源码快照、版本与已有许可证 |
| [foundations/](foundations/README.md) | 实际引用的一阶理论及 Weyl/contact 基础文件，逐字快照与 SHA-256 |
| [history/](history/README.md) | 本项目会话的用户可见消息、脱敏工具记录及导出清单 |

## 复现

```bash
git clone https://github.com/sunnowrun/gtempo-shift.git
cd gtempo-shift
python3 -m pip install -r requirements.txt
```

测试环境为 Julia 1.12.6、Python 3.11.8。`vendor/` 中保存的是普通文件，普通 clone 即可取得源码，无需初始化子模块。Julia Manifest 保留实际实验时的服务器绝对路径；首次迁移到其他路径，应在一个持久 Julia 会话中重新解析本地依赖：

```julia
include("scripts/setup.jl")
```

然后按 [最终复核复现步骤](theory_reassessment/README.md) 导出真实上游系数、生成独立参考并运行验证。首次探索的完整流程可在同一个 Julia 进程继续执行 `include("scripts/run_all.jl")`；它会覆盖同名结果，应先保留需要比较的原始数据。计算资源申请及缓存路径按自己的服务器配置；会话记录中的 SSH 别名、Slurm 作业号和 PID 是历史记录。

不需要重新运行物理模拟即可核验归档文件：

```bash
python3 scripts/verify_archive.py
python3 scripts/check_math.py
python3 theory_reassessment/scripts/check_math.py
```

完整原始会话、依赖版本和结果均保留为研究证据；早期失败控制和已被后续结论更新的解释也没有删除。原 Julia 会话交接描述见 [执行记录](theory_reassessment/docs/execution.md)，它不代表未来仍有存活的服务器任务。

## 来源与授权说明

上游 GTEMPO、QuAPI 及相关软件归各自作者所有，提交与已有许可证见 [THIRD_PARTY.md](THIRD_PARTY.md) 和 [vendor/sources.json](vendor/sources.json)。本仓库没有为第三方代码重新指定许可证。理论来源原路径及内容哈希记录在 [foundations/sources.json](foundations/sources.json)；无需访问原作者机器即可阅读已归档的九个基础文件。
