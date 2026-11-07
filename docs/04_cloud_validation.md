# 云端与工程验证：怎样积累可信证据

复杂场景与真GPU作业的具体验收关卡、能力前置、C01--C08场景阶梯、误差/性能门槛和产物合同，见 [工程测试协议](10_engineering_test_protocol.md)。对应 `tests/plans/complex_scene_matrix.json` 是计划清单，不是已执行结果。

把程序放上云端只是运行地点改变。工程验证需要同时证明数值正确、失败可识别、长期稳定、结果可复现、性能测量公平。任务 01--04 的本机数学、装配、无接触隐式推进及距离/剪切/弯曲 XPBD 已建立基线；本轮 05A 加入单元能量 XPBD，05B 的原平面接触/摩擦与完整实验仍待迁移。本轮没有新增 GPU 或远端执行。

## 当前交付与真实状态

| 项目 | 状态 |
|---|---|
| Windows MinGW Debug/Release | 05A 本轮完整 CTest 各 18/18 通过；任务 04 各 15/15 为历史记录，详情见 `09_lesson_05a.md` |
| Python 对照 | 保留前课局部/全局、48 步隐式与 32 步距离/弯曲 XPBD；05A 新增 6 个局部样例、12 次投影、32 步能量 XPBD，记录来源哈希 |
| 隐式推进与 CLI | 无接触/平面接触、最终残差、失败回滚和成功/数值失败退出已验证；任务06规则矩形场景接口与实验报告已加入 |
| XPBD 与 CLI | 分组几何/柔顺残差、预算提交与严格拒步已验证；32 步轨迹对照不是全部残差收敛的声明 |
| 05A 能量模型与结果格式 | 6 局部/12 投影/32 步对照通过；格式 2 四组诊断与乘子、两类能量分开；最后澄清配置文字后，两配置相关 2 项 CLI 检查再次通过 |
| Eigen 5.0.0 | 随工程保存离线头文件包；SHA-256 核验后解包，不发起依赖下载 |
| Linux/Windows CPU CI | 配置已写入 `.github/workflows/cpu.yml`；07C repository-only 验收脚本和 artifact 清单已加入，尚未远端运行 |
| Linux ASan/UBSan | CMake 开关与 CI job 已配置，尚未在 Linux 执行 |
| CPU Dockerfile | 已提供，尚未构建/运行镜像 |
| CUDA / GPU runner | 后续阶段，没有实现/连接/运行 |
| 长期仿真、任意拓扑场景基准与性能门槛 | 任务07B完成最小任意拓扑路径；07C已准备 CPU CI 合同；待真实云端 CPU runner 和复杂场景实现 |

## 第一阶段就可以上传的 CPU 工程

把 **`cloth_sim_cpp` 目录作为一个独立仓库的根目录**，保留隐藏的 `.github` 目录。不要把它埋在另一个仓库的子目录后期待这里的 workflow 自动生效；若选择单仓多项目，需要调整 workflow 所在位置和工作目录。

当前没有为你创建仓库、推送或购买云端资源。准备好仓库并启用 Actions 后，这个 workflow 配置了：

- Ubuntu 上 Debug、Release、ASan/UBSan；Windows 上 Debug、Release。
- 每项构建和 CTest 独立返回失败，设置超时。
- 保存 JUnit、CTest 日志和 CMake 配置日志，便于查失败。CTest 参数使用 `--output-junit ctest-results.xml`，相对测试构建目录生成，避免再次嵌套 `build/<preset>`。
- 一个单独的容器构建/运行任务。
- Linux Release 额外运行 `tools/run_ci_cpu_acceptance.py`，保存 CTest JUnit、环境元数据和 `task07_custom_mesh.scene` 报告。

CI 使用官方 [checkout](https://github.com/actions/checkout) 和 [upload-artifact](https://github.com/actions/upload-artifact) 的主版本标签，首次远端验证后再锁定审核过的完整提交 SHA。Runner 镜像和容器包也会更新；后续建立发布基线时锁定镜像 digest、依赖与工具链，不能仅凭 `ubuntu:24.04` 声称逐字节可复现。

在可使用 Docker 的机器上，从新工程根目录执行：

```sh
docker build -t cloth-cpu .
docker run --rm cloth-cpu
```

镜像构建阶段编译程序，容器启动执行 CTest，任何失败返回非零退出码。无图形窗口和 Python 运行时依赖。当前示例只演示这个路径，尚未实跑容器。

## 完整求解器之后的验收矩阵

| 层级 | 场景与指标 | 判定方式 |
|---|---|---|
| 数学 | 能量/梯度/Hessian 差分、质量、平移/转动不变性、厚度标度、J 延拓 | 小而固定的自动断言；每次变更都运行 |
| 动力学 | 静止、自由落体、固定边、受拉、强形变、接触/摩擦 | 位移/速度误差、固定点漂移、残差、穿透和耗散 |
| 鲁棒性 | 退化输入、非有限值、极端参数、线搜索失败 | 明确错误与退出码；不只检查程序未崩溃 |
| 并行一致性 | CPU 1/2/4/...线程，CPU/GPU，多次重复 | 预先规定的绝对+相对容差和物理质量指标 |
| 长跑 | 固定种子多场景、持续接触和反复拉伸 | 内存趋势、有限状态、失败数、约束误差与漂移 |
| 性能 | 无接触/接触/强形变，小中大网格 | 同模型、同精度下的总墙钟、分阶段计时与内存 |

每个测试要声明模型与离散误差预期。隐式 Euler 自带数值耗散，非凸问题不能笼统要求总能量每步必降；XPBD 距离与能量约束的残差也不能直接拿同一绝对阈值比较。

05A 的 `C=E/s` 和柔顺残差均为无量纲量，`s` 单位为 J；正柔度对应的平方势一般不同于原 NH 能量。验收需要分别核对两类能量、乘子、局部梯度和轨迹，不能仅用输出中一个叫 energy 的列判断守恒或材料等价性。05B已加入接触/摩擦的独立力学与耗散账；06增加场景哈希、提交边界和计时合同；07B已完成当前有效范围的原Python/C++交叉验收，07C将同一 C++ 合同接入 Linux runner，但远端尚未运行。 

## NVIDIA CUDA 真设备验证

先在一台明确型号的 NVIDIA GPU 上完成单卡路径；不需要一开始租多卡。云厂商与机型在真正做 CUDA 阶段时按可用性、FP64 能力与预算选择，不在当前提前承诺设备或速度。

GPU 测试可以运行在云 VM 上管理的 runner；GitHub 的 [self-hosted runner](https://docs.github.com/en/actions/concepts/runners/self-hosted-runners) 允许自选硬件和环境，需要自行维护机器。它是可选部署方式，不代表当前已经部署。

每个 GPU 作业必须记录 GPU 型号/显存、驱动、CUDA runtime/toolkit、架构目标、编译器、精度模式、源码提交、场景和种子。核启动与同步错误必须导致测试失败，不能在没有 GPU 时把 GPU 测试空跑后算通过。

[NVIDIA Compute Sanitizer](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html) 提供 memcheck、racecheck、initcheck、synccheck。后续真实可执行程序将逐项接入，并设置检查错误导致非零退出：memcheck 检查越界/对齐等，initcheck 检查未初始化的设备全局内存访问，synccheck 检查同步误用。racecheck 主要检查共享内存访问冲突，不能据其通过就证明共享顶点的全局写入无竞争；仍需用算法设计、对照测试和并发检查来验证。

当前不提供会假装 CUDA 可用的占位后端或 GPU CI 绿灯。到 GPU 阶段先实现一个可跑且可测的批量单元核，再逐步形成设备驻留的完整求解路径。

## 性能报告应保存什么

每次运行保存场景配置、种子、网格/自由度/约束数、步长、迭代上限和容差、精度、后端、线程数、硬件/系统/工具链、版本标识、是否收敛、失败原因、内存、分阶段计时及端到端墙钟。

先确定足够预热的区间，再多批次重复，报告 p50/p95 与总耗时；同时保存每个样本。GPU 同时报告设备事件时间和包含传输/同步/接触/诊断的总时间。固定记录频率和后处理设置，避免优化后只是少算了工作。

性能门槛等稳定硬件的基线形成后再设定。共享 CPU CI 主要检查正确性，不适合把小幅耗时波动当回归；GPU kernel 加速也不等于整次模拟加速。失败率和误差不合格的样本不能从报告里消失。
