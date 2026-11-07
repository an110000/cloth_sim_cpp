# 第一课：从 Python 函数走到 CMake 工程

本课只关注四件事：找到入口、构建一个库、运行程序、让测试自动判定对错。暂时不用理解 `material.cpp` 的每一行公式。

## 先运行一次

打开 PowerShell，逐行执行；上一行成功后再执行下一行：

```powershell
Set-Location 'D:\biancheng\myself\physics_animation\物理引擎\cloth_sim_cpp'
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 2
ctest --preset cpu-debug
.\build\cpu-debug\cloth_inspect.exe
```

当前机器已有 CMake 和 MinGW，所以这一步不需要安装 CUDA 或下载数学库。

第一条 CMake 命令是**配置**：读取工程描述，检查编译器，为选定构建工具生成构建规则。`-G 'MinGW Makefiles'` 指定与这台电脑 MinGW 配套的构建规则；只在首次选择生成器时需要关心它。

第二条是**构建**：调用编译器，把 `.cpp` 编译、链接成库和可执行程序。`--parallel 2` 允许同时处理两个编译任务。它不会让物理算法自动并行。

第三条是**测试**：CTest 运行登记的三个程序，任何一个返回失败退出码，整体就失败。Release 测试也使用显式检查，避免标准 `assert` 被优化配置关闭后出现假通过。

最后一条真正运行示例。应得到：

```text
vertices=9
triangles=8
total_mass_kg=1
rest_element_energy_J=0
rest_area_ratio=1
```

当前没有时间推进，所以这不是布料下落模拟。

## 只按这个顺序看四个文件

1. `apps/cloth_inspect.cpp`：C++ 的 `main()` 入口。创建材料和网格，计算质量，取第一个三角形做数学检查。
2. `include/cloth/mesh.hpp`：声明数据结构和函数。可以类比 Python 模块对外提供的接口。
3. `src/mesh.cpp`：实现网格与质量。`nx/ny` 是小格数量；顶点数为 `(nx+1)(ny+1)`，三角形数为 `2 nx ny`。
4. `CMakeLists.txt`：告诉构建系统哪些源码构成哪些目标，以及目标之间的关系。

头文件不需要“单独运行”。`main()` 调用库的函数，链接阶段把它们接在一起。

## 理解 CMake 中这几句就够了

```cmake
add_library(cloth_core src/mesh.cpp src/material.cpp)
target_include_directories(cloth_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(cloth_core PUBLIC cxx_std_17)

add_executable(cloth_inspect apps/cloth_inspect.cpp)
target_link_libraries(cloth_inspect PRIVATE cloth_core)
```

`cloth_core` 是可复用的库目标，`cloth_inspect` 是可运行的程序目标。以后测试程序、CPU 求解器和工具都可以复用同一个库。

`${CMAKE_CURRENT_SOURCE_DIR}/include` 是公开头文件的位置。`PUBLIC` 表示库本身和使用这个库的目标都需要这个目录；C++17 要求也同样传给使用者。这里公开头文件本身用了 C++ 类型，明确使用要求有助于后来加入其他目标。

`target_link_libraries(... PRIVATE cloth_core)` 声明示例程序使用物理库。CMake 因此知道先构建库，也把库的公开头文件目录传给示例程序。相比全局设置，这种写法能把将来 CUDA/测试/工具的依赖关系保持清楚。

`CMakePresets.json` 把常用配置命名为 `cpu-debug`、`cpu-release`，包括构建目录与开关。本工程使用 CMake 3.21 支持的 preset 格式。官方说明见 [CMake Presets](https://cmake.org/cmake/help/v3.21/manual/cmake-presets.7.html)。

## 把公式对应到代码

每个三角形的参考面积为 `A`，质量为 `rho * thickness * A`，平均分给三个顶点。不是把整块布料质量简单平均给所有点：边界点和内部点连接的三角形数量不同。

材料坐标 `X` 是二维的，当前位置 `x` 是三维的：

```text
Dm = [X1-X0, X2-X0]        2×2
Ds = [x1-x0, x2-x0]        3×2
F  = Ds * inverse(Dm)      3×2
C  = transpose(F) * F      2×2
J  = sqrt(det(C))          面积比
E_triangle = thickness * A * W(F)
```

`evaluate_triangle` 返回能量、对九个位置分量的梯度和 `J`。注意力是负梯度。这里尚未移植 Hessian，也没有用有限差分替代物理内核中的解析梯度；有限差分仅用于独立验证。

原 Python 的 `J<0.05` 延拓和 `det(C)<=1e-14` 退化拒绝已保留。三维中的膜 `F` 为 3×2，`J` 为非负面积比，不是三维体单元的带符号体积比。

## 你来做三个小练习

每次只做一个修改，先预测结果，再构建和运行。做完可以把预测、实际输出或疑问发到聊天里，下一次围绕你的结果讲解。

**练习 A：数顶点和三角形。** 将示例中的 `make_rectangular_mesh(2, 2, 1.0, 1.0)` 改成 `(4, 3, 1.0, 1.0)`。应输出 20 个顶点、24 个三角形，总质量仍为 1 kg。想一想：加密网格为什么不增加总质量？

**练习 B：厚度与质量。** 把示例开头材料改为 `const cloth::Material material{1.0e6, 0.3, 0.002, 0.05};`。网格不变时质量应加倍。示例用同一个材料厚度计算质量和能量，避免两处参数失配。静止能量仍为零，不能用静止能量证明刚度加倍；非零形变的厚度标度由 `test_material` 检查。

**练习 C：确认测试真的会失败。** 打开 `tests/test_mesh.cpp`，找到总质量预期为 `6.0` 的那一条，临时改成 `7.0`；重新构建后运行 CTest，应只有 `mesh_and_mass` 失败并给出原因。把预期恢复为 `6.0`，再次构建/测试，应恢复全绿。记录“修改测试源码后为什么需要重新构建”。

修改示例参数不会改变库测试里独立构造的固定场景。两类程序职责不同：示例方便观察，测试保护明确的不变量。

## Debug 与 Release

完成上述练习后恢复或保留你理解的示例参数，再体验另一个构建目录：

```powershell
cmake --preset cpu-release -G 'MinGW Makefiles'
cmake --build --preset cpu-release --parallel 2
ctest --preset cpu-release
```

Debug 方便调试，Release 开启优化。两个目录彼此独立；源码没变也需要分别构建。以后只有 Release 的同条件测量适合用来讨论性能，Debug 更适合检查问题。

如果使用 Visual Studio 多配置生成器，配置时省略 MinGW 的 `-G`，build/test preset 已指定 Debug 或 Release；程序在对应的配置子目录。不要在同一个现有构建目录混用 MinGW 和 Visual Studio。

本课完成标准：你能独立构建/测试，解释库和程序的区别，预测网格数量与质量，演示一次真实测试失败。下一课再讲解析 Hessian 和稀疏装配。
