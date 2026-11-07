# 任务 07B：任意三角拓扑输入

这一小课把规则矩形网格扩展为“场景文件直接提供三角网格”。它仍然只做布料三角膜和距离 XPBD 的最小闭环：不加入 BVH、自碰撞、CCD、复杂障碍或 CUDA。

## 1. 你要学会什么

- 用 `Mesh` 分离二维参考坐标、三维静止位置和三角形索引；
- 解释为什么参考三角形面积决定质量，而当前三维边长决定距离约束的基准；
- 用 C++17 `std::optional` 表示“使用自定义网格”与“生成规则网格”两种输入；
- 用 `std::map` 收集带编号的场景键，再检查编号从零连续，避免输入顺序影响结果；
- 用 `std::set<std::pair<size_t,size_t>>` 去重三角形的三条边，让共享边只生成一个结构约束；
- 把失败输入测试、场景 CLI 重放测试和 Release CTest 放进同一个验收报告。

## 2. 离散化合同

每个顶点有参考坐标 `X_i=(xi_i,eta_i)`、三维静止位置 `x_i` 和当前状态位置 `q_i`。对三角形 `t=(i,j,k)`，参考面积为

\[
  A_t^0=\frac12\left|(X_j-X_i)\times(X_k-X_i)\right|.
\]

质量使用 `rho h A_t^0/3` 分配给三个节点，因此孤立节点和零面积三角形必须在构造阶段拒绝。任意拓扑只改变索引集合；质量公式不依赖矩形行列编号。

对距离 XPBD，当前 07B 为每条唯一三角边生成

\[
 C_{ij}(q)=\|q_i-q_j\|-\ell_{ij},\qquad
 \ell_{ij}=\|x_i-x_j\|.
\]

这使共享边不会被重复投影。弯曲约束和非邻接自碰撞尚未加入，所以“任意拓扑已输入”不等于“任意布料物理都完成”。

## 3. 关键代码阅读

`SceneConfigInput::custom_mesh` 是可选值。场景解析器读入 `mesh_reference.N`、`mesh_position.N` 和 `mesh_triangle.N`，先检查三组编号连续，再交给 `validate_scene_config`。校验顺序是：数组长度、有限值、索引范围、三角形非退化、固定点范围。

```cpp
std::set<std::pair<std::size_t, std::size_t>> edges;
for (const auto& t : mesh.triangles) {
    for (const auto& e : {std::pair{t[0], t[1]},
                          std::pair{t[1], t[2]},
                          std::pair{t[2], t[0]}}) {
        edges.insert(std::minmax(e.first, e.second));
    }
}
```

这里的 `std::minmax` 只做无向边规范化；真正的三角形绕序和参考面积仍由网格校验负责。`std::set` 的代价是 `O(E log E)`，07B 先用它保证确定性；大型网格再测量是否换成排序后的向量或并行哈希表。

## 4. 实际验收

```powershell
Set-Location 'D:\biancheng\myself\physics_animation\物理引擎\cloth_sim_cpp'
cmake --build build/cpu-debug --parallel 4
ctest --test-dir build/cpu-debug --output-on-failure
cmake --build build/cpu-release --parallel 4
ctest --test-dir build/cpu-release --output-on-failure
python tools/run_task07b_acceptance.py
```

本轮结果：Debug 和 Release 各 29/29；原 Python `validation/` 为 71 项；07B 自定义网格场景成功提交 1 步；14 份既有参考源码散列保持一致。报告位于 `output/task07b_acceptance.json`。

## 5. 现在还不能声称什么

当前自定义输入只验证合法三角拓扑和距离 XPBD 结构边。没有实现：任意拓扑的弯曲邻接构建、缝合多裁片、OBJ/资产单位转换、BVH、连续碰撞、自碰撞、运动学障碍、刚体/流体耦合、GPU、真实云端、Linux/MSVC、sanitizer 和 Docker。下一阶段先把 07B 证据接入真实云端 CPU CI，再做 C01--C08 复杂场景；CUDA 只在 CPU 热点基线稳定后测量。

## 练习

1. 把两个三角形改成同一条边重复三次，说明为什么 `set` 仍只会产生一条距离约束，以及为什么非流形检查仍应是后续任务。
2. 把 `mesh_reference.2` 改到与 0、1 共线，观察构造阶段的错误；解释为什么不能让质量函数把这个三角形的面积当成零后继续运行。
3. 对一个 `n` 个三角形的网格估算边去重的时间复杂度，并提出一个保持确定性输出的排序实现。
