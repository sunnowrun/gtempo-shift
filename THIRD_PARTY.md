# 第三方源码与来源

`vendor/` 按实验时实际使用的提交保存六个仓库的普通文件快照，不包含这些仓库内部的 Git 数据库。准确来源见 [vendor/sources.json](vendor/sources.json)；原始实验记录见 [source_versions.txt](results/source_versions.txt)。快照保持上游源码和已有版权声明，不将其视为本项目原创代码。

| 软件 | 上游仓库 | 已有许可证文件 |
|---|---|---|
| GTEMPO | https://github.com/guochu/GTEMPO | [LICENSE](vendor/GTEMPO/LICENSE) |
| QuAPI | https://github.com/guochu/QuAPI | [LICENSE](vendor/QuAPI/LICENSE) |
| ImpurityModelBase | https://github.com/guochu/ImpurityModelBase | [LICENSE](vendor/ImpurityModelBase/LICENSE) |
| Z2Tensors | https://github.com/guochu/Z2Tensors | [LICENSE](vendor/Z2Tensors/LICENSE) |
| ExpExp | https://github.com/guochu/ExpExp | 该快照无独立 LICENSE 文件，保留原有声明 |
| FiniteMPSAlgorithms | https://github.com/guochu/FiniteMPSAlgorithms | 该快照无独立 LICENSE 文件，保留原有声明 |

Manifest 中其余 Julia 注册依赖没有复制进仓库，由 Julia 包管理器按锁定记录解析。大型运行缓存、编译产物、服务器认证配置和内部 Git 数据库不属于研究归档；研究源文件、原始结果、日志和会话导出均随仓库提供。
