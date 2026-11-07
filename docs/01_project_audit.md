# 原 Python 项目审计与迁移依据

本报告分析的原项目为 `D:/biancheng/myself/physics_animation/python/job/cloth_sim`。下文形如 `core/energy.py:45` 的位置，均相对于这个原项目；迁移文件位置则相对于新项目 `cloth_sim_cpp`。

目标是把已有的数值研究基础逐步迁移为可验证、可维护的 C++ 程序，再开展 CPU 并行、GPU 和云端测试。当前的新项目处于第一阶段：构建与局部物理内核。工程级可靠性、完整求解器和 GPU 加速仍需后续逐项建立证据。

## 1. 如何读这份报告

本报告区分三种证据，避免把读到的代码、旧结果与本轮验证混在一起。

| 标记 | 含义 | 能说明什么 |
|---|---|---|
| 本轮实跑 | 本次实际执行测试得到的结果 | 这些特定场景在当前环境下通过 |
| 静态审查 | 根据源码控制流、公式、接口作出的判断 | 发现设计和实现风险；不等于已经运行复现每个问题 |
| 历史报告 | 原项目文档中已有的实验记录 | 提供回归线索；本轮没有重新执行长时间实验 |

文件存在、能够编译、测试通过，都有各自的范围。三个 CTest 测试程序也不等于原 Python 的全部验证已经迁移。

## 2. 原项目已经具备什么

原项目的核心是一个二维材料域嵌入三维空间的三角形 Neo-Hookean 膜。它已有解析能量、梯度、Hessian、质量装配、隐式 Euler、两类 XPBD，以及地面接触和 Coulomb 摩擦。它也具备较多数学验证、数值实验和结果可视化，适合作为 C++ 迁移时的参考实现。

主要模块和数据流如下：

```text
ClothConfig
    ↓ 参数检查
create_scene → ClothMesh / 材料 / 质量 / 固定点 / 初始状态
    ↓
solver.step(x, v, fixed_nodes, dt)
    ↓
接触与摩擦处理
    ↓
能量、应力、残差、迭代次数、状态历史
    ↓
实验脚本 / 验证脚本 / 图像与动画
```

| 责任 | 原项目位置 | 迁移时需要保留的约定 |
|---|---|---|
| 参数与有效性检查 | `core/config.py:11` | SI 单位、时间步、网格尺寸、材料和接触参数 |
| 参考网格和拓扑 | `core/mesh.py:11` | `nx/ny` 计单元格，顶点数为 `(nx+1)(ny+1)`，每格两个三角形 |
| 统一场景 | `core/scene.py:36` | 初始位置在 X-Z 平面，重力沿世界 -Y；材料顶边固定 |
| 材料本构 | `core/material.py:80` | 杨氏模量、泊松比、厚度与面积塌缩延拓 |
| 单元数学 | `core/energy.py:15`、`:45`、`:54`、`:89` | 参考数据、能量、梯度、Hessian；内力为负梯度 |
| 质量装配 | `core/mass.py`、`core/mass_sparse.py:21` | 集中质量与一致质量是两种不同存储和运算路径 |
| 稠密隐式参考 | `solvers/implicit.py:20` | 小问题的参考实现，不宜用于大规模内存基准 |
| 稀疏隐式实现 | `solvers/implicit_sparse.py:43` | 批量单元计算、稀疏装配、Newton、线搜索与地面罚接触 |
| XPBD | `solvers/xpbd.py:12`、`solvers/xpbd_neohookean.py:10` | 约束逐个更新，迭代顺序影响结果 |
| 接触与摩擦 | `core/contact.py:7`、`:138` | 区分步后位置投影和求解器内部罚接触 |
| 仿真编排与记录 | `core/simulation.py:190` | 时间推进、保存频率、耗散累计、能量与应力后处理 |

一个容易混淆的细节：`ClothMesh.x` 在 `core/mesh.py:62` 使用 X-Y 平面；统一入口 `create_scene` 在 `core/scene.py:53` 重新生成 X-Z 平面初始状态。C++ 迁移应明确采用统一场景的约定，不能把两个初始位置来源混用。

### 数据布局是迁移接口的一部分

- 参考坐标 `X` 为 `N×2`，当前坐标和速度为 `N×3`，三角形索引为 `M×3`，见 `core/mesh.py:15`。
- 隐式系统按顶点交错展开自由度：`[x0,y0,z0,x1,y1,z1,…]`，固定顶点扩展为三个固定自由度，见 `solvers/implicit_sparse.py:298`。
- 原集中质量按自由度重复三次，见 `core/mass_sparse.py:64`；XPBD 转成每顶点一个质量，见 `solvers/xpbd_neohookean.py:32`。
- 新 C++ 的 `lumped_masses` 直接返回每顶点一个标量质量。比较总质量时直接求和，不能再照搬 Python 的除以 3。后续隐式装配时需要显式扩展到三个坐标。

## 3. 物理模型与求解器必须分别理解

### 3.1 单元模型

原代码采用以下运动学量，见 `core/energy.py:15`、`:33`：

```text
Dm = [X1-X0, X2-X0]       2×2
Ds = [x1-x0, x2-x0]       3×2
F  = Ds inv(Dm)           3×2
C  = FᵀF                 2×2
J  = sqrt(det(C))         面积比
```

单元弹性能为 `E_e = thickness × reference_area × W(F)`，单位为焦耳；单元质量为 `density × thickness × reference_area`，单位为千克。厚度必须在能量和质量中一致，见 `core/material.py:83`、`core/energy.py:45`。

本构使用 `mu=E/[2(1+nu)]`、`lambda=E nu/[(1+nu)(1-nu)]`，见 `core/material.py:123`；迁移时应保留这套模型，避免无意换成三维体材料的另一套 lambda 公式。

当 `J<j_min`，代码用 `ln(J)` 的二阶 Taylor 延拓；完全退化仍抛异常，见 `core/material.py:44`、`:155`。这使接近塌缩时的能量与导数保持可计算，但不提供自碰撞、拓扑修复或任意退化恢复能力。因为 `J` 来自 `FᵀF`，它也不是三维实体中带符号的体积 Jacobian。

### 3.2 三条求解路径并非仅有速度差异

| 路径 | 实际处理 | 比较时的限制 |
|---|---|---|
| 隐式 Euler + Newton | 求解惯性、膜弹性、重力；稀疏版本可把地面罚势能放入目标函数 | 稠密版本没有相同的求解内地面接触；接触后不能直接视作数值等价路径 |
| XPBD 距离/剪切/弯曲 | 多个标量约束顺序投影 | 它的距离/弯曲柔度模型与 NH 膜不相同 |
| XPBD 单元能量约束 | 每个三角形令 `C_e=E_e/energy_scale`，再进行 XPBD 投影 | 复用了 NH 能量及其梯度，但对应的柔度势能通常与原 `E_e` 不等价 |

隐式目标函数见 `solvers/implicit_sparse.py:327`。除以共同的 `dt²` 后，可以读作：

```text
Phi(x) = (x-y)ᵀ M (x-y)/(2 dt²) + E_elastic(x)
         - f_externalᵀx + E_ground_penalty(x)
y = x_old + dt v_old
```

罚接触系数为 `kappa_i=contact_eta × m_i/dt²`，并对固定点置零，见 `solvers/implicit_sparse.py:317`。摩擦仍在 Newton 之后对速度施冲量，见 `core/simulation.py:263`。因此接触罚势能与摩擦并未形成统一的联立非光滑求解。

**静态审查与公式推导：XPBD-NH 的模型差异。** `core/constraints.py:106`、`:127` 明确使用 `C=E/scale`，并以 `alpha_tilde=compliance/dt²` 更新乘子。对于正柔度 `alpha`，标准 XPBD 标量柔度势能为：

```text
U = C²/(2 alpha)
  = E²/(2 alpha scale²)
```

这通常不等于原始 NH 势能 `E`，其回复力也不同。若 `alpha=0`，则是硬约束极限，不能将上式直接除以零。把 `energy_compliance` 调小，是更强地要求该单元回到零能量状态，不能直接解释为与原材料一致的杨氏模量。

`core/simulation.py:49` 对 XPBD 状态计算同一 NH 能量用于图表比较。这是共同的形变度量；对距离 XPBD，甚至对当前定义的能量 XPBD，都不能据此宣称图中的 NH 能量就是其柔度模型的实际势能。后续若要做“同一物理模型、不同求解器”的实验，需要先明确约束参数化和能量等价关系，再讨论速度与误差。

## 4. 本轮实际验证

**本轮实跑：46 项通过，用时 16.83 秒。** 执行了以下核心集合，未运行全量绘图、长时间收敛研究和大规模动画：

```powershell
python -m pytest -q -p no:cacheprovider `
  validation/test_core_pytest.py `
  validation/test_contact.py `
  validation/test_collapse.py `
  validation/test_sparse_equivalence.py `
  validation/test_gravity_drop.py `
  validation/test_xpbd_neohookean.py
```

这些测试覆盖能量/梯度/Hessian 的基本一致性、质量与场景约定、接触基本行为、塌缩延拓、稠密与稀疏装配/轨迹等价、重力下落、XPBD-NH 基本行为。它们没有覆盖全部参数、全部接触边界或硬件环境。

执行环境：Python 3.9.10，NumPy 1.25.2，SciPy 1.10.1，pytest 6.2.5。本机还能找到 CMake 4.0.0-rc5 和 MinGW GCC 14.2.0；PATH 中未找到 Ninja、MSVC `cl`、Clang 或 CUDA `nvcc`。这只是当前工具发现结果，不代表已经确认机器没有 NVIDIA 显卡。

原报告中关于能量峰值、空间/时间收敛、接触修复前后对比的数值属于**历史报告**。例如 `验证分析报告.md:218` 已说明，非凸隐式 Euler 的总能量不增并无一般证明；特定场景的能量结果不能推广为所有输入的保证。

## 5. 迁移前应记入任务清单的风险

### 5.1 测试失败不一定使自动流程失败

**静态审查：** `validation/test_contact_in_newton.py:29` 使用全局失败计数，`:32` 的 `check()` 只累计和打印，不抛异常；直到脚本入口 `:119` 才按失败数退出。通过 pytest 收集这些 `test_*` 时，`check(False)` 不会自动使该用例失败。后续应改成直接断言或统一异常，并验证“故意引入错误时 CI 确实变红”。本轮 46 项测试未包含这个文件，不能以这 46 项作为该接触回归的运行结果。

**静态审查：** `validation/test_energy_conservation.py:384`、`validation/test_accuracy_cost.py:188` 捕获仿真异常后打印并继续。能量分析在不足两个成功方法时正常 `return`，见 `validation/test_energy_conservation.py:393`。这些研究脚本尚不能直接作为失败即终止的云端验收程序。

梯度、Hessian 和空间收敛等文件还含独立 `main()`、绘图及交互显示，例如 `validation/test_gradient.py:107`、`validation/test_hessian.py:128`。应拆成无图数值测试与可选报告生成，避免将“文件名叫 test”误认为“已经被自动测试执行”。

### 5.2 未收敛的状态会继续被接收

**静态审查：** 稀疏 Newton 在线搜索无法接受候选步时 `break`，最终仍返回当前状态与 `converged=False`，见 `solvers/implicit_sparse.py:448`、`:451`、`:465`。`core/simulation.py:247` 随即接收新位置和速度，没有针对失败状态执行拒步、回滚或减小时间步。

记录层的 `converged` 还默认取真，见 `core/simulation.py:186`；`save_every>1` 时只保存取样步的求解状态，中间未收敛步可能不会出现在保存结果中，见 `core/simulation.py:301`。工程迁移应逐步累计失败次数，并明确“接受、重试、终止”的策略。

**静态审查：** `solvers/xpbd.py:155` 和 `solvers/xpbd_neohookean.py:62` 将“位置和速度有限”作为 `converged`。这与“约束残差低于容差”含义不同，甚至可能同时出现有限状态、退化跳过和无穷残差。C++ 中应分开报告 `finite_state`、`residual`、`tolerance_satisfied`、`degenerate_count`，避免一个布尔值承担多种含义。

### 5.3 固定点位于地面下方时的投影冲突

**静态审查，尚未另行运行复现：** `core/contact.py:63` 至 `:69` 排除了固定点参与冲量计算，但 `:70` 的位置投影对所有穿透顶点执行。若固定点初始就在地面下方，它仍可能被移动到地面，违反固定位置约束。现有正常悬挂场景不会触发这一特定条件。

迁移时应先定义固定约束与障碍物冲突的语义：拒绝冲突场景，或允许固定点穿透而保持固定；再用测试落实该决定，不能依赖默认场景恰好没有冲突。

### 5.4 计时范围和比较对象尚不统一

**静态审查：** `benchmark.py:24` 只记录单次运行中各步的平均值，没有预热、重复批次、分位数、内存或失败计数；`:40` 默认仅五步。`:51` 稀疏版带地面罚接触，`:53` 稠密版不带。前几步未触地时比较尚有共同区间，延长到触地后便不是同一问题。

**静态审查：** `core/simulation.py:247` 至 `:254` 的 `step_times` 只计 `solver.step`，不包含之后的摩擦、接触投影、应力后处理和保存；`:307` 仅保存当前取样步的时间，不是保存区间的耗时之和。把这些值求和不能得到完整仿真时间。

**历史报告：** `验证分析报告.md:158` 承认部分墙钟时间受同时运行的 CPU 任务干扰。后续应提供明确命名的求解时间和端到端时间，固定场景与误差要求，记录暖机、重复次数、线程数、编译配置、硬件和收敛失败数，再报告加速比。

## 6. 如何映射到新的 C++ 项目

迁移应围绕数学责任与数据接口组织，逐项对齐原实现。

| Python 来源 | C++ 迁移位置/方向 | 验收重点 |
|---|---|---|
| `core/mesh.py` + `core/scene.py` 的网格约定 | 已有 `include/cloth/mesh.hpp`、`src/mesh.cpp` | 拓扑、面积、X-Z 初始坐标；未来场景层再加入固定边和速度 |
| `core/mass_sparse.py` 的集中质量 | 已有 `lumped_masses` | 每顶点标量质量、正质量、总质量 `rho t A` |
| `core/material.py` + `core/energy.py` 的能量/梯度 | 已有 `include/cloth/material.hpp`、`src/material.cpp` | 单位、厚度标度、刚体变换、有限差分、延拓区域 |
| Python 小单元参考 | 已有 `tools/export_python_reference.py`、`tests/data/python_reference.txt` | 可检查的输入输出、源文件摘要、普通/延拓/旋转/不同参考三角形 |
| 单元 Hessian 与全局稀疏装配 | 后续任务 | 对称性、梯度差分、装配模式、Python 对齐 |
| 隐式 Euler 与 Newton | 后续独立求解器模块 | 固定点、线搜索、退化拒绝、未收敛策略、轨迹对齐 |
| 接触与摩擦 | 后续独立接触模块 | 穿透、摩擦锥、耗散、固定约束优先级 |
| XPBD 约束族 | 后续可选求解器模块 | 先决定模型语义，再验证投影算法与残差 |
| `core/simulation.py` 与实验输出 | 后续无界面 CLI 与结果记录 | 每步状态、失败退出、参数与环境元数据、端到端计时 |
| NumPy 批量计算和顺序约束循环 | 后续 CPU 并行/GPU 后端 | 数据竞争、求和误差、传输成本、真实设备执行证据 |

CPU 并行可以先针对互相独立的单元局部计算。全局梯度/Hessian 装配会让相邻单元同时写同一顶点，不能直接给循环加并行指令。XPBD 在 `solvers/xpbd.py:129` 逐约束修改位置，改成并行 Jacobi 或图着色也会改变迭代方式，需要单独的收敛与精度验证。

GPU 规划同样需要先确定数据驻留、装配或 matrix-free 运算、线性求解策略。只在 GPU 上计算单元能量，而每步往返传输并在 CPU 求解，不足以证明整体仿真加速。

## 7. 第一阶段能力边界与学习顺序

第一阶段代码已提供网格、集中质量、NH 单元能量/梯度、异常输入检查，以及 CMake/CTest 的构建和测试入口。它是后续求解器可以调用的局部物理库；`cloth_inspect` 用于检查这些基本量。

当前尚无 C++ 时间推进、全局 Hessian、稀疏 Newton、接触摩擦、弯曲、自碰撞、OpenMP、CUDA 或完整仿真性能基准。第一阶段测试只能证明当前已实现部分的性质，不能证明布料仿真已经完整迁移。本报告的“46 项通过”特指原 Python 项目的本轮测试，C++ 构建结果另见项目交付说明。

对已学过基础 C++、刚接触 CMake 的学习者，建议先能解释以下三件事，再进入下一阶段：

1. `cloth_core` 是库，`cloth_inspect` 和测试程序是调用它的可执行文件；CMake 为这些目标描述源文件、头文件路径、语言标准和链接关系。
2. “配置、编译、测试”是三个动作：配置生成构建规则，编译产生库和程序，CTest 执行已经注册的测试。配置成功不等于编译或测试成功。
3. 一个三角形输入参考坐标和当前位置，输出能量与梯度；梯度是能量对坐标的导数，内力为它的负值。先用这个小问题建立正确性证据，再增加全局装配和时间推进。

下一阶段应以单元 Hessian、全局装配和无接触的串行隐式求解为主要学习任务。只有串行基线的模型、误差与失败策略清楚之后，CPU/GPU 加速与云端大规模测试的结果才有可解释性。
