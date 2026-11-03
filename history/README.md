# 项目会话历史

这里保存本项目当前会话的用户可见历史，包括最初问题、第一阶段答复、要求回到 GTEMPO 理论基础的修正、理论复核过程、最终结论及公开仓库请求。

- [conversation.txt](conversation.txt)：按时间排列的可读文本。
- [conversation.jsonl](conversation.jsonl)：每行一条用户或助手消息，含时间、角色和回复阶段。
- [tools.jsonl](tools.jsonl)：工具调用及结果，保留执行上下文并作凭据脱敏。
- [export_manifest.json](export_manifest.json)：原会话文件 SHA-256、快照截止时间、消息数及导出规则。
- [export_conversation.py](../scripts/export_conversation.py)：可复核的导出程序。

这是截至发布前导出时的会话快照；GitHub 创建、推送及最终发布确认可能发生在快照截止时间之后。仅导出这个项目会话，没有收集其他无关会话。内部指令、隐藏推理、压缩上下文和遥测未纳入公开历史；读取会话内部文件的工具记录保留时间和名称，内容省略。认证凭据按导出规则脱敏。数值方法、推导、源代码及执行日志另有完整归档。

会话中的本地路径、旧作业编号、暂时性错误和早期结论是当时记录，不是当前操作指令；请以最终理论复核报告与磁盘脚本为准。
