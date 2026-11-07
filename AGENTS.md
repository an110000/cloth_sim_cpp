# 项目学习约定

GPU/云端与复杂场景的工程验收协议见 `docs/10_engineering_test_protocol.md`，计划矩阵 `tests/plans/complex_scene_matrix.json`。C01--C08均为planned，不能当已执行测试。阈值先经应用需求与参考研究标定，再测性能；要求GPU的作业无设备必须环境失败，不得静默CPU回退。复杂资产、BVH、自碰撞/CCD仍在原项目迁移完成后的相应扩展阶段。CPU真实云端CI可在07先执行，不能把配置存在记为云端通过。

本项目按用户要求分任务迁移，不一口气完成所有后端。用户会基础 C++，CMake 从头学习；GPU 目标为 NVIDIA CUDA。

当前任务 01--07C 已实现当前有效范围（05A 能量 XPBD、05B 原平面接触/摩擦、06 场景/实验合同、07 迁移验收、07B 任意三角拓扑输入、07C repository-only CPU CI 合同），本机 MinGW Debug/Release 完整 CTest 各 29/29 通过；原 Python validation 71 项通过。远端 Linux runner 尚未执行；下一步是上传 artifact 后审查，再做复杂场景和 CUDA。07B 的规则与自定义网格场景、唯一结构边和失败校验已落盘；07C 的 JUnit/环境/场景报告脚本已落盘。C01--C08、GPU/云端仍未实跑。

用户最新要求：不加入 Rust；先完成原布料项目的 C++ 重写与验收，再扩展刚体、流体、通用碰撞、摩擦和多物理场耦合。总路线以 `docs/05_multiphysics_architecture.md` 为准，`02_learning_roadmap.md` 保留布料数值细节。05A/05B/06/07/07B分别保留独立讲义与源码快照。实现到07B覆盖当前有效范围，不等于复杂资产、CPU并行、CUDA和云端已经通过。

XPBD 的有限状态、按组柔顺残差和是否提交必须分开；默认预算模式可提交有效近似，严格模式残差不达标回滚。32 步原 XPBD 轨迹对照不意味着每步残差收敛。原 Python 退化时可返回零角/零梯度，C++ 明确标无效并拒步，属于有记录的保护性差异。任务 04 讲义入口 `../cloth_sim_lectures/task04.tex`，运行指南 `docs/08_lesson_04.md`，快照 `../cloth_sim_lectures/code/task04/`；不覆盖历史任务快照。

05A 复用 `Constraint` variant 与 `XpbdSolver`，增加能量约束和第四组残差；NH 诊断能量与平方约束势分开。05B 在隐式路径加入步内地面罚势和步后摩擦，在 XPBD 路径加入步后平面响应及残差复测。固定点不投影、分离速度保留、法向冲量记录实际速度改变；非零 beta/restitution 仅投影路径支持。数值手册第 3 版同步推导模型、罚刚度、动能账与接触分裂的性质；总册第 7 版，不改旧任务快照。

Eigen 5.0.0 头文件包由本机已有源码离线重打包，来源和许可证见 `third_party/README.md`；CMake 使用 `file(SHA256)` 核验后 `file(ARCHIVE_EXTRACT)`。不要把本地包冒称已经独立核验的官方发布归档，不在正常配置中下载依赖。Linux、MSVC、sanitizer、Docker、GPU 尚未实跑。

每次独立的迁移或学习任务都要同步在相邻 `../cloth_sim_lectures/` 中提供或更新分层 LaTeX 学习讲义，包含：理论推导、离散化、真实关键代码片段、相关 C++ 特性、CMake、执行与验收、练习及答案提示。任务涉及 GPU/云端时，补充真实设备操作、错误检查、计时边界与证据。

当前 C++ 求解器保持无窗口 headless；`tools/render_experiment.py` 是独立分析层，把已提交的 `CLOTH_EXPERIMENT_TEXT` 帧渲染为 PNG/GIF 和诊断图。它不得重新推进物理或把图片通过当成数值验证；图片必须和文本报告、场景哈希一起保存。

复用总册和每任务单册结构。保存对应源码快照与来源，不因工程更新静默改写旧任务快照。明确区分已实现、实际验证、预期操作和计划。更新讲义后编译并检查排版；不要把没有真实 GPU 的空跑计为通过。

日常小修正更新所属任务讲义即可，不必每次创建新任务编号。继续使用 `cloth_sim_lectures/templates/task_template.tex` 的章节组织方式。

独立数值分析学习手册位于 `../cloth_sim_numerics/`，根入口 `numerical_analysis.tex`。用户要求整个开发过程持续推导数值算法性质；每新增或实质改变算法，同步更新该手册的假设、误差/一致性/收敛/稳定/条件数/守恒耗散、实际证据与局限。明确线性模型定理不自动适用于非凸膜、接触和流体；数值预算提交不等于收敛。更新实现后显式核对本册源码基线，不自动改写历史分析或旧任务源码快照。
