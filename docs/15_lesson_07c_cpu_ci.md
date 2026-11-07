# 任务 07C：真实云端 CPU CI

07B 证明了自定义三角网格在本机可运行。07C 把这条验收路径接入 Linux runner：构建、CTest、任意拓扑场景、环境元数据和 JUnit 必须一起保存。当前仓库已经写入 workflow，但只有远端 runner 真正上传报告后，才能登记为云端通过。

## 1. CI 和本机测试的区别

本机 Debug/Release 29/29 是开发机证据；云端 CPU CI 还要记录：

- checkout 后的提交身份和工作区状态；
- Linux、CPU 架构、CMake、编译器和 Python 版本；
- CTest JUnit、CTest 最后日志和配置日志；
- 07B 自定义网格场景的输入 SHA-256、顶点/三角形数量和提交步数；
- CTest 或场景失败时的完整产物。

`tools/run_ci_cpu_acceptance.py` 不依赖原 Python 工程，因此可以在独立 C++ 仓库中运行。它只负责仓库内 CPU 验收；原 Python 71 项交叉对照仍由本机 `run_task07b_acceptance.py` 完成。

## 2. Workflow 结构

`.github/workflows/cpu.yml` 保留 Debug、Release、sanitizer 和 Windows 矩阵。Linux Release 在普通 CTest 后再运行：

```yaml
- name: Repository CPU acceptance
  if: matrix.os == 'ubuntu-24.04' && matrix.preset == 'cpu-release'
  run: |
    python3 tools/run_ci_cpu_acceptance.py \
      --build "${{ github.workspace }}/build/${{ matrix.preset }}" \
      --output "${{ github.workspace }}/build/${{ matrix.preset }}/ci-acceptance"
```

`if: always()` 的 artifact 步骤必须放在测试之后。这样失败的 Linux、Windows、sanitizer 和容器任务也会留下 JUnit、最后日志和配置记录，不会只剩一个红色状态图标。

## 3. 接受脚本的判断

脚本执行 `ctest --output-on-failure --output-junit`，解析 JUnit 的 `tests/failures/skipped`，再运行 `tests/scenes/task07_custom_mesh.scene`。场景报告必须满足：状态为 `committed`、请求步数等于提交步数、第一帧有 4 个顶点、场景文件含 4 个节点和 2 个三角形。它还会保存场景哈希和 CTest JUnit 哈希。

```powershell
python tools/run_ci_cpu_acceptance.py `
  --build build/cpu-release `
  --output output/task07c-ci-local
```

本机可以用这条命令模拟报告生成；它不是 Linux 云端通过，因为本机的 `environment` 和 `source_identity` 不能替代远端 runner。远端报告必须从 workflow artifact 下载后再审查。

## 4. 失败分类和证据边界

构建失败、CTest 失败、场景数值失败、环境不可用和 artifact 上传失败要分别看。Docker 镜像构建通过只说明容器中的编译/CTest路径通过；它不证明云 runner、GPU 或复杂场景通过。没有原 Python 仓库时，CI 报告明确写 `repository-only`，不能把缺失的跨项目对照改写成绿色。

云端 CPU 也不等于 Linux 数值定理。它只把相同 C++ 离散问题放到另一套工具链和内核上运行；跨平台位置差异要用预先冻结的绝对/相对容差和物理指标解释。计时字段在本任务只作为环境观察，不作为重放判据。

## 5. 练习

1. 将自定义场景中的一个三角形索引改成 4，观察场景报告和 workflow 的退出码；说明为什么这属于输入失败而不是 CTest 通过。
2. 删除 `ci-ctest-results.xml` 后重新运行脚本，解释为什么即使自定义场景通过，整体报告也应失败。
3. 比较本机和云端报告中的 `compiler`、`machine`、`junit_sha256`，指出哪些字段能用于身份追溯，哪些字段不能推出性能结论。
