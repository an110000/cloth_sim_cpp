# 任务 07：原项目迁移验收

任务 07 把原 Python 验证集合、C++ Release CTest、三种任务 06 场景和参考来源散列合并成一份验收报告。

运行：

    python tools/run_task07_acceptance.py

报告写入 `output/task07_acceptance.json`，当前验收包含：

- 原 Python `validation/`：71 项通过；
- C++ Release CTest：28/28 通过；
- 距离 XPBD、隐式 Euler、能量 XPBD 三种场景入口；
- 14 份 Python 参考源码 SHA-256 逐项一致。

报告仍明确列出任意拓扑、C01--C08、GPU、云端、Linux、MSVC、sanitizer 和 Docker 未执行。当前结果是迁移验收证据，不是复杂物体或真实设备验收。

下一阶段可以把该脚本接到真实云端 CPU CI；云端报告必须保留源码、场景、构建器、JUnit、Python 日志、失败样本和环境信息。GPU 作业需要单独的真实设备与 Compute Sanitizer 证据。
