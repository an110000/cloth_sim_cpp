# 任务 06：场景与实验输出

本课把已经通过局部测试的求解器变成可复现的场景运行入口。核心新增两层：`SceneConfig` 负责构造后不可变的场景值，`ExperimentRecord` 负责提交边界、失败前缀和原子实验报告。下一项是任务 07 原项目迁移验收。

## 运行一个场景文件

场景文件采用不依赖第三方库的 `key=value` 文本格式，注释以 `#` 开头。最小例子：

```text
base=teaching
solver=xpbd-distance
mass_mode=lumped
dt=0.01
steps=8
ground_enabled=true
ground_y=-0.03
friction=0.3
static_friction=0.5
```

构建后运行：

```powershell
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 1
ctest --preset cpu-debug
.\build\cpu-debug\cloth_scene_run.exe --scene tests\scenes\task06_teaching.scene --output tmp\task06-run.txt
```

报告是固定顺序的 UTF-8 文本，包含模型、质量模式、步长、请求步数、已提交步数、计时和逐帧状态。重复运行同一场景时，去除墙钟与求解计时字段后，物理状态和报告结构应逐字节一致；计时本身受系统调度影响。输出目标已存在时程序拒绝覆盖，避免把一次不同运行伪装成同一个结果。

## 提交边界和失败

每个物理步先由求解器返回候选状态，只有 `info.committed=true` 才调用 `append_committed_frame`。因此失败步之前的帧仍是有效前缀，失败候选不会进入轨迹。完成所有请求步后才能标记 `committed`；数值失败写 `failed_numerical` 并退出2；参数/文件错误退出1。后续云端还会使用 `failed_memory`、`failed_environment`、`timeout`、`interrupted` 和 `unsupported`，不能把它们都折叠成普通数值失败。

`RunMetadata` 中的源码、配置、场景、设备和容器哈希目前使用显式 `<unset>` 占位，不能把占位符当成已经发现的环境信息。真实发布脚本必须在运行前填入哈希、编译器、后端、设备、容器和随机状态。`SceneConfig` 会拒绝非法尺寸、有限性错误、重复固定点、错误 SHA-256 格式，以及隐式接触中不适用的恢复系数/稳定化参数。

## 学习顺序和任务边界

先读 `include/cloth/scene.hpp`，理解 value type、构造校验和只读 getter；再读 `include/cloth/experiment.hpp`，沿着 `append -> finish -> validate -> atomic rename` 追踪生命周期。CMake 把两个核心源文件和两个单元测试加入 `cloth_core`/CTest；CLI 只承担文本解析和求解器编排，不把文件格式污染到物理内核。

本课已经实现配置化的规则矩形教学场景、初始位置/速度覆盖、三种求解器选择和固定平面参数。任意三角拓扑、复杂资产、BVH、自碰撞、CCD、CPU 并行、CUDA 和云端仍不属于本课；这些能力需要先完成任务 07 的原项目验收，再按路线扩展。
