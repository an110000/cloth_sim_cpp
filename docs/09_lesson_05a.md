# 第 05A 课：单元能量 XPBD

本轮迁移原 `C = E / energy_scale` 路径，继续复用串行 `XpbdSolver`。任务 05 尚未全部完成，下一项 05B 是原平面接触与摩擦，再完成 06--07 的场景和迁移验收。

工程讲义：[TeX 入口](../../cloth_sim_lectures/task05a.tex)，输出 `../../cloth_sim_lectures/output/pdf/task05a_energy_xpbd.pdf`；性质分析见 `../../cloth_sim_numerics/numerical_analysis.tex` 第 9 章。代码快照在 `../../cloth_sim_lectures/code/task05a/`，历史快照不刷新。

新增 `NeoHookeanEnergyConstraint` 保存节点、二维参考坐标、按值材料、柔度、能量尺度和乘子；`Constraint` variant 增加第三种类型，组别追加 `neohookean_energy`，原三个组索引保持不变。直接调用已有解析能量/梯度并除以尺度，不计算 Hessian。极小正尺度使用逐分量除法，避免先求倒数溢出后把合法零梯度变成 NaN。

令尺度 s 的单位为 J，C 无量纲，能量柔度 alpha 单位 1/J。正柔顺模型的实际势为各单元 `E_e² / (2 alpha_e s_e²)` 的和，一般不等于原 NH 膜能量。参考态附近可出现四次能量与零切线刚度；固定每单元参数细化网格也会改变该模型。原定义保留，未悄悄替换成另一材料。

从工程目录运行：

```powershell
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 1
ctest --preset cpu-debug
cmake --preset cpu-release -G 'MinGW Makefiles'
cmake --build --preset cpu-release --parallel 1
ctest --preset cpu-release
./build/cpu-release/cloth_xpbd_simulate.exe --model energy --output output/task05a
./build/cpu-release/cloth_xpbd_simulate.exe --model energy --energy-scale 3 --energy-compliance 1.1111111111111111e-9 --output output/task05a_scaled
./build/cpu-release/cloth_xpbd_simulate.exe --model energy --strict --iterations 1 --energy-tolerance 1e-14 --output output/task05a_strict
```

默认距离路径保留。能量参数只在 `--model energy` 下接受；能量路径不额外加入距离/剪切/弯曲约束。默认 energy scale=1 J、compliance=1e-8 /J、扫描30次，仍是固定教学场景；完整原场景不是本轮范围。Python 原能量求解器默认扫描10次，逐项对照显式使用相同扫描预算，不声称两个默认配置完全相同。

格式2输出 `config.json`、`trajectory.csv`、`diagnostics.csv`、`multipliers.csv`、`summary.json`。诊断表四组分别使用 m、m、rad、无量纲；NH 构型诊断能量与平方势分列，后者在距离模型或 hard alpha=0 时为 null。仅提交步写入轨迹与乘子；失败 summary 可描述被拒绝候选，返回的状态仍为原状态。

本机 MinGW Debug/Release **各18/18 CTest通过**，新增6局部/12投影/32步 Python 对照；导出器自身最大梯度差分相对误差2.72e-10，参考轨迹最大柔顺残差4.24e-9。数据再次生成字节一致；修改临时副本的期望值+0.01被读取器拒绝。完整测试日志在 `output/task05a_validation/`；最终CLI说明文字澄清后，两种配置相关CLI测试另行通过。

Release教学场景8步全提交，max|R|=2.6093e-7，小于默认阈值1e-6；成对缩放位置差最大9.605e-14 m、速度差最大8.399e-12 m/s。严格1次扫描、1e-14容差第一步拒绝，退出码2；参数/IO错误为1，完成请求步数为0。退化几何拒步是有记录的保护性差异，原Python可能跳过后继续。

数值手册运行 `python tools/analyze_task05a.py` 后，独立重算R与报告差最大1.50e-16，并观察到位置形式驻点失配 `s_x=x-y-WJᵀlambda` 最大约0.103 m。s_x是方程失配，不能解读为真实轨迹误差；严格模式目前仅检验柔顺R，未实现完整非线性驻点停止准则。

下一课学习原平面接触、Coulomb摩擦、冲量与耗散记账。本轮没有新增GPU、云端或其他平台实跑。
