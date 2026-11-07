# 任务 05B：原平面接触与摩擦

2026-11-04 完成三条布料路径的固定 +Y 平面接触。下一项是 **06 场景与实验输出**，之后 07 原项目迁移总体验收；继续使用 C++ / CMake / 后续 CUDA，不加入 Rust。

学习入口：`../../cloth_sim_lectures/task05b.tex`，独立 PDF 为 `../../cloth_sim_lectures/output/pdf/task05b_contact_friction.pdf`。数值分析第三版的新章 `../../cloth_sim_numerics/chapters/11_contact_friction.tex` 推导误差、稳定性、条件数与耗散边界。历史任务快照不覆盖，本课冻结到 `cloth_sim_lectures/code/task05b/`。

## 三条更新路径

| 路径 | 法向处理 | 切向处理 | 收敛判据 |
| --- | --- | --- | --- |
| 隐式 Euler，两种质量矩阵 | Newton 目标中加入罚势 `Uc=.5 sum(kappa pen²)`，`kappa=eta*m_lumped/dt²`；残差及 Hessian 同步加入 | 位置收敛后以 `Jn=dt*kappa*pen` 作节点 Coulomb 修正；没有第二次法向处理或位置投影 | 罚接触位置方程，摩擦后速度不被声称满足一个联立摩擦 Newton 系统 |
| 距离/剪切/弯曲 XPBD | 材料投影后仅投影自由节点，按实际 `m*delta_vn` 记法向冲量 | 同一节点 Coulomb 映射 | 重算**接触后**分组柔顺残差；严格模式不达标整步回滚 |
| 能量 XPBD | 同上；仍不含额外弯曲模型 | 同上 | 同上；不把 NH 诊断能量误当实际平方约束势 |

隐式法 `max_penetration_before/after` 是最终位置上摩擦阶段两端（相等），罚法允许小穿透；XPBD 是位置投影两端。固定点穿透单独报告，既不移动固定点，也不伪装成自由接触。隐式 `normal_kinetic_change=0` 仅表示没有独立法向后处理，不表示 Newton 过程没有法向运动。

## 审查后有意修正的原行为

1. 原地面函数会投影固定点；新版固定位置保持原值，穿透作为冲突诊断。
2. 原函数将向上离地速度清零；新版保留分离速度，只施加非负必要冲量。
3. 原法向冲量包含未施加到速度的 beta 项；新版冲量等于实际动量变化。beta 默认从旧 `.2` 改为 `0`；显式正 beta 可增加法向动能，单独记账。
4. 原 `np.isclose` 还含隐式相对容差；新版使用明确的绝对 proximity tolerance。
5. 原 consistent 路径使用 `diag(M)` 作为接触/摩擦质量；新版使用物理行和节点质量。P1 三角形当前 `m_lumped=2*diag(M)`，因此旧 consistent 参考轨迹采用 **mu=0、eta_old=30 → eta_port=15**。这是参数映射对照，不能称同参数轨迹回归；新 consistent 摩擦另有独立测试。
6. 零切向速度的旧两个摩擦入口有不同计数；新版正冲量且零切速统一计粘着。原审计不要求这些旧计数逐值相同。

节点摩擦耗散 `lumped_friction_dissipation_J` 用物理集中质量计算；`actual_friction_kinetic_change_J` 使用求解器真实质量矩阵。consistent 的非对角交叉项使两者不同，真实动能甚至可能增加。当前为有记录的节点分裂近似，尚非联立广义质量摩擦接触算法。单元测试有真实三节点布料质量增能反例。

## 构建与运行

```powershell
cmake --preset cpu-debug -G 'MinGW Makefiles'
cmake --build --preset cpu-debug --parallel 1
ctest --preset cpu-debug
cmake --preset cpu-release -G 'MinGW Makefiles'
cmake --build --preset cpu-release --parallel 1
ctest --preset cpu-release
./build/cpu-release/cloth_simulate.exe --ground-y -.03 --contact-eta 30 --friction .3 --static-friction .5 --output output/contact_implicit
./build/cpu-release/cloth_xpbd_simulate.exe --ground-y -.03 --friction .3 --static-friction .5 --output output/contact_distance
./build/cpu-release/cloth_xpbd_simulate.exe --model energy --ground-y -.03 --friction .3 --static-friction .5 --output output/contact_energy
```

XPBD 可显式设置 `--restitution` 和 `--contact-beta`；隐式接触拒绝这些不适用项。接触参数须搭配 `--ground-y`；静摩擦须不小于动摩擦，所有参数必须有限。省略地面继续旧的无接触场景。这里只扩充固定六节点教学 CLI，通用配置入口留在 06。

`contact.csv` 只写已提交时间步，单位见表头；XPBD 另写 `contact_split.csv`。summary 保留最后一次尝试的候选指标；严格失败时不能把它们当作已提交轨迹或累计耗散。已有 diagnostics.csv 的列顺序保持；新增配置/summary 字段为追加。NH 诊断在距离模型中若不适用或溢出，写 null 并标记 `diagnostic_energies_valid=false`，不让附加诊断中断已提交的距离模型状态；该标记进入后续验收门槛，有限预算提交本身不代表质量合格。

## 本机证据

- Windows、MinGW GCC 14.2.0、CMake 4.0.0-rc5，Debug/Release 完整 CTest 各 **23/23**。审查增加巨大有限速度零摩擦与非适用诊断两项断言后，Debug 相关 **2/2** 再次通过；最终 Release 完整验收包含它们。
- 新增公共内核、隐式接触、XPBD 接触、原参考读取器、CLI 合同五项测试；覆盖解析解、罚势/完整目标差分、冲量/摩擦锥/耗散、分离速度、固定穿透、后接触残差、几何退化和回滚。
- 原参考：4 个摩擦、5 个地面、2 个罚势、3 项旧行为审计；四轨迹各 8 步，共 **32 步真实触地**。重导出 fixture 和来源清单逐字节一致，14 份原 Python 来源散列不变。原隐式最终残差最大 `7.90e-11`；局部罚梯度差分误差最大 `9.45e-10`。
- 独立真实 Release CPU 脚本 `cloth_sim_numerics/tools/analyze_task05b.py`：能量 XPBD 8 步接触前 R 最大 `2.56097e-7`、接触后 `7.68437e-6`，显示一次分裂会扰动材料约束；两配置严格模式已覆盖最终拒步。距离路径前/后分别 `8.36619e-7 m` / `2.33733e-4 m`。
- 一致质量 CLI 同一摩擦阶段累计实际动能变化 `-0.00468447 J`，节点耗散 `0.00925918 J`，不能混写。另有独立两节点 SPD 教学反例与引擎三节点反例，不把教学矩阵数据当引擎测量。
- 固定 k 的单自由度活动罚振子教学实验最后观测阶 `0.99943`；它只验证同一个光滑活动分支 ODE，不是连续布料或含冲击接触的收敛证明。固定 eta 时模型随 dt 变化的实验单列。

物理来源仍是原 `core/contact.py`、`core/friction.py`、`solvers/implicit_sparse.py` 和 `core/simulation.py`，有意修正见上。法向/切向冲量的拓展阅读为 [Catto GDC 2006 原讲义](https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf)；目前没有迁移其中的通用刚体多接触求解器。

Linux、MSVC、sanitizer、Docker、CUDA、复杂碰撞与云端尚未实跑。测试协议矩阵所有复杂场景仍为 planned；已有脚本/CI 配置不会被登记成云端通过。
