# 任务 02：解析 Hessian、两种质量与全局装配

本任务已完成。当前可以构造全局弹性能、梯度、CSR Hessian、质量矩阵和固定自由度子系统；还没有时间推进或 Newton 求解。下一项是任务 03：隐式 Euler / Newton。

完整推导与真实代码片段见相邻 `../cloth_sim_lectures/task02.tex`，PDF 为 `../cloth_sim_lectures/output/pdf/task02_hessian_assembly.pdf`。

## 先运行这个小例子

在 C++ 工程根目录逐行执行：

```powershell
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 2
ctest --preset cpu-debug
.\build\cpu-debug\cloth_assembly_inspect.exe
```

示例是一平方米的单格网格，两个三角形、四个节点。共有 12 个位置自由度；固定节点 2、3 后，剩下 6 个自由自由度。预期完整 Hessian 为 `12x12`，约化矩阵为 `6x6`，总质量约 1 kg。

`hessian_stored_entries=126` 指存储项数，包含显式零。两个局部 9×9 矩阵共生成 162 项；共享节点 0、3 的 36 个重复位置合并后，剩 126 项。一致质量存储 42 项；`mass_00_kg=1/6`、`mass_03_kg=1/24`，节点质量应取行和而不是只取对角线。

## 本次代码入口

| 文件 | 要理解的责任 |
|---|---|
| `material.hpp/.cpp` | `evaluate_triangle_with_hessian`，解析方向导数；原一阶接口继续保留 |
| `sparse.hpp/.cpp` | COO 重复项合并、CSR 不变量、查询与矩阵乘法；不是线性求解器 |
| `mass.hpp/.cpp` | 集中/一致质量，统一输出交错坐标的 3N×3N 矩阵 |
| `assembly.hpp/.cpp` | 只装配膜弹性，局部索引映射和共享贡献累加 |
| `dofs.hpp/.cpp` | 自由/固定映射、子向量、A_FF 与零固定增量提升 |
| `apps/cloth_assembly_inspect.cpp` | 最小可观察示例 |

`lift_from_free` 只用于增量，不是绝对位置重建；非零 Dirichlet 绝对值问题仍需单独修正右侧。全固定情况是合法 0×0 子系统；重复固定点明确拒绝。

## 验收与边界

2026-11-04，本机 MinGW GCC 14.2.0 / CMake 4.0.0-rc5，Debug、Release 各 7/7 CTest 通过。检查包含解析不变量、单元/全局差分、两种质量、固定映射、非法输入，以及 11 个 Python 单元 Hessian、3 个全局系统。

将实验副本中的一个参考 Hessian 条目加 1，测试返回 1 并报告差异；正式参考数据保持不变。原 Python 源码没有修改。新参考生成器记录原模块和旧输入数据哈希，正常 C++ 构建与 CTest 不依赖 Python。

暂未实跑 Linux、MSVC、sanitizer、Docker 或 GPU。当前装配每次建立/排序 COO，没有优化成固定模式原位更新；稠密导出只供小型对照。Hessian 对称不代表正定，没有完整模拟性能或加速倍数结论。

## 三个小练习

1. 将示例固定顶点 `{2,3}` 改成 `{0}`，预测约化矩阵大小，再构建运行：应为 9×9。
2. 将一致质量改成集中质量：存储项变成 12，`M00=1/3`、`M03=0`，总质量和行和不变。
3. 用自己的话解释“能量相加、梯度散加、Hessian 三元组合并”，并说明为何 `H_FF` 现在不能直接被称为完整 Newton 矩阵。

每次只改一处，先写预测，再运行并解释结果。实验结束后恢复示例，保留记录。下一讲把质量、惯性、外力与 Newton 残差接到这个系统上。
