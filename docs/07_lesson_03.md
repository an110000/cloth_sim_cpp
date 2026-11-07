# 任务 03：无接触隐式 Euler、Newton 与失败回滚

任务 03 已实现第一条串行时间推进路线。固定场景的无窗口示例、质量与重力、稀疏 LU、下降方向检查、Armijo 回溯与最终残差共同构成一个可验收的时间步；原 XPBD、弯曲、地面接触和摩擦仍待后续迁移。

完整分层推导见相邻 `../cloth_sim_lectures/task03.tex`；从连续动量方程消去速度，得到

```
y = x_n + dt*v_n
Phi(x) = 0.5*(x-y)^T*M*(x-y) + dt^2*(E(x)-f_g^T*x)
r = M*(x-y) + dt^2*(grad(E)-f_g)
A = M + dt^2*H
A_FF*d_F = -r_F
```

`Phi` 的单位是 kg·m²，`r` 是 kg·m，`A` 是 kg。材料厚度与时间步是不同参数。质量采用 `Material::thickness`，避免弹性和惯性各用一份不一致的厚度。

## 先跑固定的小场景

在 C++ 工程根目录，使用现有 MinGW 构建目录：

```powershell
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 1
ctest --preset cpu-debug
.\build\cpu-debug\cloth_simulate.exe --output output/task03_lumped
.\build\cpu-debug\cloth_simulate.exe --mass consistent --output output/task03_consistent
```

默认是 6 个节点、4 个三角形，上边三个节点固定，下边预先轻微拉伸；E=20000 Pa、nu=0.3、厚度 0.001 m、密度 800 kg/m³，dt=0.01 s，共 8 步。两种质量都应打印 `steps_completed=8/8 converged=true failure_reason=none`。`--parallel 1` 控制编译，并非仿真并行。

每个输出目录包含：

- `config.json`：完整初始位置/速度、拓扑、材料、固定节点和求解设置；当前 CLI 尚不从此文件加载任意场景。
- `trajectory.csv`：初始帧及每个已提交帧的位置/速度，共 9×6=54 条数据。
- `summary.json`：完成步数、有限状态、收敛、提交、失败原因和最后尝试的残差。

失败示例：

```powershell
.\build\cpu-debug\cloth_simulate.exe --max-iterations 0 --output output/task03_failure
$LASTEXITCODE
```

应返回 **2**，完成 0 步，原因 `iteration_limit`；CSV 只有 6 个初始节点。参数/IO 错误返回 **1**，正常完成返回 **0**。`finite_state=true` 仍可能是回滚后的旧状态，必须同时检查收敛和提交。失败报告里的残差描述最后被评估的候选点，并不一定属于已回滚的输出状态。没有自动减半时间步重试。

## 读代码时只抓四件事

1. `SimulationState` 只保存位置和速度；`ImplicitEulerSolver` 拥有参考几何、材料、质量和常重力。输入是 const 引用，返回独立结果。
2. 重力使用完整 `f_g=M*g`。预测 `y+dt²*g` 只是初值，不能在收敛结果上再补一次重力。
3. Eigen SparseLU 处理一般方阵；不要求矩阵正定，也不保证 Newton 方向下降，所以仍检查 `r·d<0` 和 Armijo。位移只改变方向，不改变被验收的方程。
4. 最后一次允许的更新后仍重算残差。只有最终残差合格且输出有限才提交；失败返回完整原状态。

`iterations` 统计尝试的 Newton 更新，`accepted_steps` 统计被接受的更新；两者不同。固定位置取本步输入位置，成功时固定速度置零。输入固定速度若非零，一致质量惯性耦合仍保留，详见讲义。

## 离线依赖与 CMake

`third_party/eigen-5.0.0-headers.tar.gz` 来自本机已有独立 Eigen 5.0.0，保留上游许可证。它是本地重打包，不是声称经过独立核验的官方发行包；精确内容以 SHA256 锁定，见 `third_party/README.md`。CMake 校验并在构建目录解压，无需网络，也无需原电脑路径。

`cloth_eigen_headers` 是 INTERFACE 头文件目标，`cloth_core` PRIVATE 使用它；公开物理头文件不暴露 Eigen 类型。当前只构建串行 Eigen 路径。每次方向尝试仍转换、分析和分解矩阵，尚未复用符号分解或拓扑数值槽位。

## 实际验证

2026-11-04，Windows / MinGW GCC 14.2.0 / CMake 4.0.0-rc5：Debug、Release 各 11/11 CTest 通过。新测试覆盖：

- 一般不定稀疏系统、奇异失败与空系统。
- 静止、两种质量的手算自由落体、固定位置/速度、目标差分与独立最终残差。
- 最后一次允许更新后的收敛、压缩膜实际使用位移方向、塌缩预测回退。
- 预算耗尽、线搜索失败、溢出与无效输入；失败回滚。
- 原 Python 6 条轨迹共 48 步对照；生成器另对其中 8 步做原稠密/稀疏对照。
- CLI 两种质量成功、预算失败、JSON/CSV 与退出码合同。

改坏实验副本的期望位置会使参考测试返回失败；正式基线不变。正常 C++ 构建/测试不依赖 Python。Linux、MSVC、sanitizer、Docker 与 GPU 尚未实跑，云端配置不能当成已通过的结果。

## 这一课的三个练习

1. 自由落体 h=0.1、g=-10、v0=0，一步位移为何是 -0.1 而不是 -0.05？从后向 Euler 的位置式解释时间离散误差。
2. 将同一输出目录实验改为 `--mass consistent`，先预测固定边是否移动、轨迹是否必须相同，再运行。固定位置不变，变形轨迹一般不同。
3. 运行零预算失败例，解释为什么 `finite_state=true` 与 `committed=false` 可以同时成立；确认 CSV 没有伪造出 t=0.01 的帧。

下一课任务 04：迁移原距离/剪切 XPBD 与弯曲。新刚体、流体和通用碰撞仍等待整个 C++ 迁移验收。
