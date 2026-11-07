# Eigen header dependency

`eigen-5.0.0-headers.tar.gz` is a **local repackaging**, not the upstream release archive.
It contains the unmodified `Eigen/` header directory and upstream license/readme files
copied on 2026-11-04 from the existing local checkout
`D:/biancheng/myself/physics_animation/femdemo/external/eigen`.
`Eigen/Version` reports `EIGEN_VERSION_STRING "5.0.0"` with no prerelease suffix.
No upstream commit identity has been independently verified; the exact packaged
bytes are pinned by SHA256:

`b5e679981cf1c843ade160aeda435e05de917f653b29e8d5ab854acdfa128adf`

The build needs only this repository, C++17 and CMake 3.21+. It neither downloads
dependencies nor reads the original local source directory. CMake verifies this
archive and extracts it under each build directory. The `EIGEN_MPL2_ONLY` build
definition restricts the Eigen modules selected by our code; upstream license
files remain in the archive. `EIGEN_DONT_PARALLELIZE` keeps this lesson serial.

We currently use Eigen's square sparse LU with pivoting rather than an SPD-only
solver, because a membrane Newton matrix may be indefinite. This first adapter
rebuilds/factorizes each attempt; it does not yet reuse symbolic analysis or claim
optimal memory/performance. Changing the archive requires review, a new hash and
rerunning numerical/reference tests.

Upstream: https://gitlab.com/libeigen/eigen

Sparse solver documentation: https://libeigen.gitlab.io/eigen/docs-5.0/group__TopicSparseSystems.html
