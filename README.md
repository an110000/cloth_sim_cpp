# cloth_sim_cpp：多物理场学习项目的 C++ 数值基础

这是配合学习逐步完成的迁移项目。原 Python `cloth_sim` 保留为数学参考；新目录独立构建。当前 **任务 01--07C 已实现当前有效范围**：05A 为单元能量 XPBD，05B 为三条路径的原平面接触/摩擦，06 为配置化场景与实验输出，07 为原项目迁移验收，07B 增加合法任意三角拓扑输入和距离 XPBD 结构边，07C 增加 repository-only Linux CPU CI 验收脚本和 workflow。远端 runner 尚未执行；复杂场景和 CUDA 仍未完成。

05B 本机 MinGW Debug/Release 完整 CTest 各 **23/23** 通过；06 合同加入后三种求解器路径各 28/28；07B 后 Debug/Release 各 **29/29**，并实际运行原 Python 71 项验证、07B 自定义网格场景和 14 项来源散列。consistent 对照使用显式质量/罚刚度映射，不称同参数一致；GPU、云端与复杂场景未新增实跑。

采用 C++ / CMake 与后续 CUDA。**先完整重写现有布料项目并验收，再扩展刚体、流体、通用碰撞、摩擦和双向耦合。** C++ 统一承担物理状态、几何与数值求解；CUDA 加速测量后的批量热点。本目录保留现有名称与可运行基线，不把计划模块登记为已实现。

总设计和14个分步任务见[迁移优先的总路线](docs/05_multiphysics_architecture.md)，配套分层TeX讲义在相邻`../cloth_sim_lectures`。当前完成到07C的有效范围；下一步是远端 CPU runner 实跑与 artifact 审查，再做复杂场景和 CUDA，再扩展新物理模块。

## 从这里开始

配套新增 [数值算法性质分析手册](../cloth_sim_numerics/README.md)（独立 TeX/PDF），从当前 FEM、隐式 Euler、Newton 和 XPBD 推导误差、收敛、稳定与耗散；后续每个算法同步补充性质和实验，明确条件结论与实际观测的区别。

1. 阅读 [项目分析](docs/01_project_audit.md)，弄清楚原项目已有和缺少什么。
2. 浏览 [扩展总路线](docs/05_multiphysics_architecture.md)；[布料数值子路线](docs/02_learning_roadmap.md) 保留原有 FEM/GPU 学习细节。
3. 按顺序跟做 [第一课](docs/03_lesson_01.md)、[第二课](docs/06_lesson_02.md)、[第三课](docs/07_lesson_03.md)、[第四课](docs/08_lesson_04.md)、[05A](docs/09_lesson_05a.md)、[05B：平面接触与摩擦](docs/11_lesson_05b.md)、[06：场景与实验输出](docs/12_lesson_06.md) 和 [07B：任意拓扑](docs/14_lesson_07b_arbitrary_topology.md)。每课保留练习和独立讲义；数值分析手册同步推导模型、误差与收敛边界。
4. [云端与工程验证](docs/04_cloud_validation.md) 区分已配置、已本机验证和待真实设备验证。

## 当前实现

| 模块 | 内容 | 验证 |
|---|---|---|
| `cloth_core` | 规则三角网格、集总质量 | 拓扑、质量守恒、输入检查 |
| `evaluate_triangle` | 3×2 Neo-Hookean 膜的能量、解析梯度、面积比 | 有限差分、刚体不变性、厚度缩放 |
| `evaluate_triangle_with_hessian` | 普通/延拓区解析 9×9 Hessian | 梯度差分、原始对称性、平移零模、转动协变 |
| `SparseMatrix` / 质量 | COO→CSR、集中/一致质量、矩阵乘法 | 重复合并、行和守恒、标度、异常输入 |
| `assemble_elastic` / `DofMap` | 全局能量/梯度/Hessian、固定自由度约化 | 全局差分、小网格对照、合法全固定空系统 |
| `solve_sparse_lu` | Eigen 5.0.0 稀疏 LU 与后向误差检查 | 非正定/奇异系统、非法输入与数值失败 |
| `ImplicitEulerSolver` | 隐式 Euler、Newton、可选步内地面罚势与步后摩擦、失败回滚 | 两种质量、原轨迹、固定点、最终残差、接触与真实质量矩阵动能变化 |
| `Constraint` / `XpbdSolver` | 距离、双对角剪切、解析二面角梯度；顺序 Gauss-Seidel 投影 | 手算、差分、原约束顺序与 32 步轨迹；有限状态、分组残差和提交分别报告 |
| `NeoHookeanEnergyConstraint` | `C=E/energy_scale`，复用 NH 解析能量/梯度和串行 XPBD | 05A 新增 6 个局部样例、12 次投影和 32 步原能量路径对照；尺度补偿、固定点、退化与失败语义见第 05A 课 |
| `contact.hpp` | +Y 平面接触、Coulomb 摩擦、罚势及逐节点冲量/耗散 | 原有效域对照与旧行为审计；固定点不投影，分离速度保留，法向冲量对应实际速度改变 |
| Python 对照 | 11 个原实现生成的固定样例 | 普通形变、J 延拓及连接点、任意参考形状、反向绕序 |
| `cloth_inspect` | 无窗口的最小检查程序 | 输出网格数、质量和静止能量 |
| `cloth_assembly_inspect` | 单格四节点的装配检查 | 12×12 完整、6×6 约化系统，质量矩阵和存储项 |
| `cloth_simulate` | 固定教学场景的无窗口推进与结果输出 | 成功退出、数值失败非零退出与状态回滚；不是任务 06 全部场景接口 |
| `cloth_xpbd_simulate` | 固定场景的两种 XPBD 模型、可选平面接触、预算/严格模式 | 四组诊断、乘子与两类能量；接触前后残差重新报告，接触可能破坏已投影约束 |
| `SceneConfig` / `cloth_scene_run` | 版本化 key=value 场景、不可变校验、三种求解器编排 | Debug/Release 场景合同、重复物理帧、失败前缀和原子实验报告 |
| 工程配置 | CMake/CTest、Debug/Release、Linux sanitizer、CPU CI、CPU Dockerfile | 本机 MinGW Debug/Release 已通过；其他环境待运行 |
| 任意拓扑场景 | `SceneConfigInput::custom_mesh`、严格网格校验、唯一三角边距离约束 | Debug/Release 29/29、自定义网格场景；弯曲、自碰撞、复杂障碍待做 |
| CPU CI 验收 | Linux workflow、CTest JUnit、环境元数据、自定义网格报告、失败 artifact | 本机脚本可模拟；远端 runner 尚未执行 |
| 可视化适配器 | 将 C++ 实验文本报告渲染为 PNG、GIF 和速度/位移诊断图 | 已用任务 06 三帧场景生成并检查图片；不改变 headless 求解器 |

单位统一为 SI：位置 m、质量 kg、杨氏模量 Pa、厚度 m、能量 J。接口中的 `gradient` 是 `dE/dx`，力为它的负值。`positions` 使用原场景的 X-Z 平面，默认重力沿 -Y；原 `core/mesh.py` 单独构造时的 X-Y 嵌入与场景约定不同，这里采用场景约定。

## 在当前电脑运行

需求：C++17 编译器、CMake ≥ 3.21、相应构建工具。当前电脑已有 MinGW，PowerShell：

```powershell
Set-Location 'D:\biancheng\myself\physics_animation\物理引擎\cloth_sim_cpp'
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 2
ctest --preset cpu-debug
.\build\cpu-debug\cloth_inspect.exe
```

预期输出：9 个顶点、8 个三角形、总质量 1 kg、静止单元能量 0、面积比 1。`--parallel 2` 只并行编译，不表示仿真已经多线程。

任务 02 新示例：`.\build\cpu-debug\cloth_assembly_inspect.exe`。完整讲义入口为相邻目录的 `task02.tex`，PDF 在其 `output/pdf/task02_hessian_assembly.pdf`。

任务 03：`.\build\cpu-debug\cloth_simulate.exe --steps 8 --mass lumped --output output/task03`。见 [第三课](docs/07_lesson_03.md)；其独立讲义保留无接触版本。当前 05B 可选接触见新课，固定教学场景仍不代表任务 06 已迁移。

任务 04：`.\build\cpu-debug\cloth_xpbd_simulate.exe --output output/task04`。见 [第四课](docs/08_lesson_04.md)；TeX 入口 `../cloth_sim_lectures/task04.tex`，PDF 为其 `output/pdf/task04_xpbd_bending.pdf`。默认预算模式可提交有效近似，不能据 `finite_state` 或已提交就称残差收敛；`--strict` 可拒绝最终分组柔顺残差不达标的步。退化约束在 C++ 中明确标为无效并拒步，区别于原 Python 返回零角/零梯度的行为。

任务 05A：`.\build\cpu-debug\cloth_xpbd_simulate.exe --model energy --output output/task05a`。见 [05A 指南](docs/09_lesson_05a.md)；独立入口 `../cloth_sim_lectures/task05a.tex` 保留该课版本。`--energy-scale` 单位 J、`--energy-compliance` 单位 J 的倒数，`--energy-tolerance` 无量纲。输出中的 NH 诊断量 `E` 与实际平方约束势 `E²/(2*compliance*energy_scale²)` 分开。

任务 05B：`.\build\cpu-debug\cloth_xpbd_simulate.exe --model distance --ground-y -.03 --friction .3 --static-friction .5 --output output/task05b`。隐式 CLI 使用 `--ground-y` 和 `--contact-eta` 启用步内罚接触；详见 [本轮指南](docs/11_lesson_05b.md) 与 `../cloth_sim_lectures/task05b.tex`。计划矩阵的 C01--C08 全部仍为 planned，当前 +Y 平面能力不代表复杂物体碰撞已实现。

任务 06：`.\build\cpu-release\cloth_scene_run.exe --scene tests\scenes\task06_teaching.scene --output output\task06-run.txt`。场景文件使用严格的 `key=value` 文本，不依赖 JSON 库；输出报告把已提交帧和失败前缀分开，计时字段独立记录。详见 [任务06指南](docs/12_lesson_06.md) 与 `../cloth_sim_lectures/task06.tex`。

任务 07B：`.\build\cpu-release\cloth_scene_run.exe --scene tests\scenes\task07_custom_mesh.scene --output output\task07b-run.txt`。`mesh_reference.N`、`mesh_position.N` 和 `mesh_triangle.N` 必须从零连续编号；参考三角形必须有正面积。当前距离 XPBD 对每条唯一三角边生成一个结构约束，报告和限制见 [07B 指南](docs/14_lesson_07b_arbitrary_topology.md)。

任务 07C：`python tools/run_ci_cpu_acceptance.py --build build/cpu-release --output output/task07c-ci-local`。这会生成仓库内 CPU 验收报告；它不替代远端 Linux runner。workflow 位于 `.github/workflows/cpu.yml`，说明见 [07C 指南](docs/15_lesson_07c_cpu_ci.md)。

可视化验证：`python tools/render_experiment.py --report output/visual/task06_teaching_report.txt --scene tests/scenes/task06_teaching.scene --output output/visual/task06_teaching`。输出 `frame_*.png`、`animation.gif` 和 `diagnostics.png`；它读取 C++ 已提交帧，不在 Python 中重新推进物理。当前 C++ 核心仍为无窗口程序，图片和动画属于分析/验收层。

Eigen 5.0.0 头文件包随工程保存，配置时先核验 SHA-256，再解包到构建目录；不联网，也不依赖本机原来的 Eigen 路径。来源和许可见 [依赖说明](third_party/README.md)。

Linux 上从工程目录执行：

```sh
cmake --preset cpu-release
cmake --build --preset cpu-release --parallel 2
ctest --preset cpu-release
./build/cpu-release/cloth_inspect
```

Windows 使用 Visual Studio 时，在它的开发者终端中配置，省略 MinGW 的 `-G`；多配置生成器的程序一般在 `build/cpu-debug/Debug/`。已经配置过的构建目录不要切换生成器；切换时使用新的构建目录。

## 目录

```text
include/cloth/      网格、本构、CSR、装配、状态、隐式与 XPBD 推进接口
src/               CPU 串行实现
apps/              最小无窗口程序
tests/             自动失败的 C++ 测试与 Python 固定对照数据
tools/             从原 Python 工程导出参考值
docs/              分析、路线、任务01--05、测试协议与云端验证说明
third_party/       离线 Eigen 头文件包、许可证与来源哈希说明
.github/workflows/ 云端 CPU 验证配置，尚未在远端执行
```

历史记录保留：任务 04 为 Debug/Release 各 15/15，05A 为各 18/18；其参考轨迹和独立单册不改写。本轮 05B 另增加 32 个实际触地步，其中 consistent 对照以 `eta_old=30→eta_port=15、mu=0` 显式映射，局部摩擦/罚势与旧行为审计另列。有限状态、残差达标与提交分开；本轮数值结果见 [05B 指南](docs/11_lesson_05b.md)。没有新增 Linux、MSVC、sanitizer、Docker、GPU 或云端实跑，也不据此宣称完整重写验收。

新工程的正常编译、CTest 不需要 Python。只有主动更新参考数据时才需要原项目和 NumPy；更新后必须审查差异，不能用“重生成期望值”消除失败：

```powershell
python tools/export_python_reference.py --source 'D:\biancheng\myself\physics_animation\python\job\cloth_sim'
```

参考数据保存原文件 SHA-256 与 NumPy 版本，方便追溯。
