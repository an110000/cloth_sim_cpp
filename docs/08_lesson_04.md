# 任务 04：距离、剪切 XPBD 与二面角弯曲

这一课迁移原 `solvers/xpbd.py` 的无接触距离路径：先预测位置，按结构→剪切→弯曲的原顺序投影，再重建速度。对应分层讲义在相邻 `../cloth_sim_lectures/task04.tex`，PDF 为 `output/pdf/task04_xpbd_bending.pdf`。旧三课的源码快照不覆盖。

## 先做一个能解释的实验

在 C++ 工程根目录逐行执行：

```powershell
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 1
ctest --preset cpu-debug
.\build\cpu-debug\cloth_xpbd_simulate.exe --output output/task04
```

场景是 2×1 格、6 节点、4 三角形，上边三个节点固定。默认每步 0.01 秒、8 步、30 轮扫描；结构 7 条、剪切 4 条、弯曲 3 条。原二维参考坐标决定距离静止长度，三维静止位置决定弯曲静止角；动态初态另有轻微变形。

输出目录包含完整 `config.json`、初始及已提交位置/速度 `trajectory.csv`、分组指标历史 `diagnostics.csv`、最终 `summary.json`。正常完成有 54 条节点数据，另有表头。当前 CLI 保存固定场景配置，不读取任意配置文件；完整原场景/实验迁移仍属于任务 06。

## 这一课最重要的区别

对于正柔度，势能为 `U=C*C/(2*compliance)`。距离约束 `C` 的单位是 m，角度约束为 rad；两类柔度的单位不同。一步内的标量更新为：

```
alpha_tilde = compliance / dt^2
delta_lambda = (-C - alpha_tilde*lambda) /
               (sum(inverse_mass_i * |gradient_i|^2) + alpha_tilde)
delta_x_i = inverse_mass_i * gradient_i * delta_lambda
```

每个物理时间步重置乘子，步内各轮累计。固定节点必须用零逆质量，原实验脚本通过 `mass[fixed]=inf` 达到同样效果；C++ 保留有限物理质量，在 `step()` 内处理固定集合。仅投影后把固定点改回去，会错误地保留它在分母中的可移动质量。

原始几何偏差 `C` 和柔度指标 `R=C+compliance/dt²*lambda` 不同。软约束可能有非零 C 而 R 很小；R 只是逐约束指标，不是完整全局 Newton 或速度约束验收。不能混合米与弧度取一个统一的最大误差。

## 亲自对比预算与严格模式

```powershell
.\build\cpu-debug\cloth_xpbd_simulate.exe --iterations 1 --output output/task04_budget_one
.\build\cpu-debug\cloth_xpbd_simulate.exe --strict --iterations 1 --distance-tolerance 1e-14 --angle-tolerance 1e-14 --output output/task04_strict
$LASTEXITCODE
```

本机第一条完成 8 步并返回 0，但最终 `projection_converged=false`，如实表示有限预算近似。第二条返回 2、完成 0 步，原因 `projection_not_converged`，轨迹只有 6 条初始节点数据。两者的返回状态都可以有限。

`finite_state`、`projection_converged`、`committed` 分别表达有限数组、分组柔度指标与提交政策。默认预算模式允许有效有限的近似，严格模式要求指标达标。严格失败时状态回滚，报告指标与乘子描述最后候选；程序没有自动重试。参数/IO 错误返回 1，完成所选政策返回 0，数值或严格拒步返回 2。

## 代码入口与本课 C++

| 文件 | 阅读责任 |
|---|---|
| `state.hpp` | 公共位置/速度数据；隐式与 XPBD 接口复用 |
| `constraints.hpp/.cpp` | `std::variant` 两类约束、`std::visit`、解析角度梯度与局部投影 |
| `xpbd.hpp/.cpp` | 步内顺序更新、乘子重置、固定权重、分组指标与回滚 |
| `cloth_xpbd_simulate.cpp` | 固定场景配置、已提交轨迹、诊断、退出码 |

弯曲构建用 map 找到边，用独立 vector 保留首次出现顺序；不能直接按 map 键排序输出。CMake 继续使用 C++17 和同一个 `cloth_core`，新增两个实现、CLI 与四组测试；没有下载新依赖或实现多线程/GPU。

## 已验证与有意修正

2026-11-04，Windows / MinGW GCC 14.2.0 / CMake 4.0.0-rc5：Debug、Release 各 **15/15 CTest** 通过。新证据包括手算软/硬投影、二面角解析梯度差分、不变量、拓扑与几何检查、局部溢出回滚、自由落体、乘子重置、提交政策、CLI 文件与退出码。

原 Python 对照包含 7 个局部案例、14 次连续投影、2 组有序构建、4 条轨迹共 32 步，比较位置/速度/乘子和分组指标。32 步表示轨迹移植对照，不表示所有步都达到默认阈值。三组参考最大柔度指标约为：结构 8.3662e-7 m、剪切 1.8904e-7 m、弯曲 4.3991e-4 rad。

参考数据再生逐字节一致；修改实验副本中一个期望值会真实返回失败，正式基线不变。普通 C++ 构建与测试不依赖 Python。

原代码在退化几何可能返回零角/零梯度，C++ 明确标记无效并拒步；原构建器略过非双面边，C++ 拒绝非流形连接。有效参考案例保留原顺序与算法，保护性差异单独记录。距离网络不自动等价于 NH 膜材料。Linux、MSVC、sanitizer、Docker、云端和 GPU 尚未实跑。

## 三个小练习

1. 用 `--no-bending` 重跑，预测计数和轨迹：结构/剪切仍为 7/4，弯曲变 0，轨迹可变。
2. 同样一轮预算，在两种政策下检查完成步数、物理时间与 CSV：未达严格阈值不能写新帧。
3. 解释为何每轮清零乘子会改变算法，以及软约束的非零几何偏差为何不自动说明失败。

下一课任务 05：能量 XPBD 与原地面接触、摩擦。完成 C++ 重写与验收后，再扩展刚体、流体和通用碰撞。
