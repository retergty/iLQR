# Matrix 与 Eigen 替换兼容方案

## 目标

本方案用于将 `Tests` 之外的项目源码逐步从 Eigen 迁移到 `Matrix` 目录下的自定义固定尺寸矩阵库，同时保留测试代码与 Eigen 的兼容能力。

迁移目标分为三层：

1. 核心求解器、模型、动力学、积分、控制器等非 `Tests` 源码默认使用自定义矩阵类型。
2. `Tests` 目录暂不强制改写，仍可使用 `Eigen::Matrix`、`Eigen::Dynamic`、LU、SVD 等测试辅助能力。
3. 通过兼容层允许测试中的 Eigen 固定尺寸对象与核心源码中的自定义矩阵对象互相转换，降低一次性替换风险。

## 当前依赖现状

项目大多数核心类型已经通过 `Types.hpp` 暴露：

```cpp
template <typename Scalar, int Rows, int Cols>
using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

template <typename Scalar, int Rows>
using Vector = Eigen::Matrix<Scalar, Rows, 1>;
```

这是迁移的主要切入点。只要非测试源码都继续依赖 `Matrix` 和 `Vector` 别名，而不是直接依赖 `Eigen::Matrix`，就可以通过替换 `Types.hpp` 完成大部分后端切换。

但测试工具中存在较多 Eigen 专属能力，例如：

- `Eigen::Dynamic`
- `resize()`
- `fullPivLu()`
- `lu().solve()`
- `Eigen::LDLT`
- `Eigen::JacobiSVD`
- 块矩阵拼接与动态 KKT 求解

因此短期不建议移除测试侧 Eigen。

## 推荐文件结构

建议新增或调整如下文件：

```text
Types.hpp
Matrix/MatrixEigenLikeApi.hpp
Matrix/EigenCompat.hpp
```

职责划分：

- `Types.hpp`：项目公共线性代数类型入口，负责选择 Eigen 后端或自定义矩阵后端。
- `Matrix/MatrixEigenLikeApi.hpp`：为自定义矩阵补充项目需要的 Eigen 风格固定尺寸接口。
- `Matrix/EigenCompat.hpp`：在测试或过渡期启用 Eigen 互转构造、赋值与转换。

## `Types.hpp` 后端选择建议

`Types.hpp` 应继续提供全局 `Matrix` 和 `Vector` 名称，但不再直接固定到 Eigen：

```cpp
#pragma once

#if defined(ILQR_LINALG_BACKEND_EIGEN)

#include <Eigen/Core>

template <typename Scalar, int Rows, int Cols>
using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

template <typename Scalar, int Rows>
using Vector = Eigen::Matrix<Scalar, Rows, 1>;

#else

#include "Matrix/Matrix.hpp"
#include "Matrix/Vector.hpp"
#include "Matrix/SquareMatrix.hpp"

#if defined(ILQR_ENABLE_EIGEN_COMPAT)
#include <Eigen/Core>
#endif

namespace ilqr_linalg {

template <int Dim>
struct StaticDim {
  static_assert(Dim >= 0, "Custom Matrix only supports fixed dimensions.");
  static constexpr std::size_t value = static_cast<std::size_t>(Dim);
};

template <typename Scalar, int Rows, int Cols>
struct MatrixSelector {
  using type =
      matrix::Matrix<Scalar, StaticDim<Rows>::value, StaticDim<Cols>::value>;
};

template <typename Scalar, int Rows>
struct VectorSelector {
  using type = matrix::Vector<Scalar, StaticDim<Rows>::value>;
};

#if defined(ILQR_ENABLE_EIGEN_COMPAT)
template <typename Scalar>
struct MatrixSelector<Scalar, Eigen::Dynamic, Eigen::Dynamic> {
  using type = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
};

template <typename Scalar, int Cols>
struct MatrixSelector<Scalar, Eigen::Dynamic, Cols> {
  using type = Eigen::Matrix<Scalar, Eigen::Dynamic, Cols>;
};

template <typename Scalar, int Rows>
struct VectorSelector<Scalar, Eigen::Dynamic> {
  using type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
};
#endif

}  // namespace ilqr_linalg

template <typename Scalar, int Rows, int Cols>
using Matrix = typename ilqr_linalg::MatrixSelector<Scalar, Rows, Cols>::type;

template <typename Scalar, int Rows>
using Vector = typename ilqr_linalg::VectorSelector<Scalar, Rows>::type;

#endif
```

这样可以保证：

- 正式源码默认使用自定义固定尺寸矩阵。
- 测试侧如需 `Eigen::Dynamic`，可在启用兼容宏时继续落到 Eigen 动态矩阵。
- 若需要快速回退，可通过 `ILQR_LINALG_BACKEND_EIGEN` 临时恢复 Eigen 后端。

## 自定义矩阵需要补充的最小接口

为了减少非 `Tests` 源码改动，自定义矩阵建议补齐一组 Eigen-like API。

### 矩阵接口

```cpp
static Matrix Zero();
static Matrix Ones();
static Matrix Constant(Type value);
static Matrix Identity();

constexpr static size_t rows();
constexpr static size_t cols();
constexpr static size_t size();

bool isApprox(const Matrix& rhs, Type eps = Type(1e-8)) const;
bool isZero(Type eps = Type(1e-8)) const;
void setConstant(Type value);

template <size_t R, size_t C>
auto topLeftCorner();

template <size_t R>
auto topRows();

template <size_t R>
auto bottomRows();
```

这些接口可以直接复用已有的 `setZero()`、`setIdentity()`、`setAll()`、`slice<P, Q>()` 等实现。

### 向量接口

```cpp
Type squaredNorm() const;

template <size_t N>
auto head();

template <size_t N>
auto tail();
```

`squaredNorm()` 可直接委托给已有的 `norm_squared()`。

## Eigen 兼容层建议

在 `ILQR_ENABLE_EIGEN_COMPAT` 打开时，固定尺寸自定义矩阵支持从 Eigen 固定尺寸对象构造与赋值。

建议形式：

```cpp
#if defined(ILQR_ENABLE_EIGEN_COMPAT)
template <typename Derived>
Matrix(const Eigen::MatrixBase<Derived>& other) {
  static_assert(Derived::RowsAtCompileTime == M ||
                Derived::RowsAtCompileTime == Eigen::Dynamic);
  static_assert(Derived::ColsAtCompileTime == N ||
                Derived::ColsAtCompileTime == Eigen::Dynamic);

  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      (*this)(i, j) = other(static_cast<int>(i), static_cast<int>(j));
    }
  }
}

template <typename Derived>
Matrix& operator=(const Eigen::MatrixBase<Derived>& other);

operator Eigen::Matrix<Type, M, N>() const;
#endif
```

注意事项：

- 兼容层只建议用于测试、过渡或验证，不应成为核心求解路径的常态依赖。
- 生产构建默认不定义 `ILQR_ENABLE_EIGEN_COMPAT`，避免核心代码重新依赖 Eigen。
- 动态尺寸仍交给 Eigen，不建议在当前阶段给自定义矩阵库加入动态尺寸能力。

## 非 `Tests` 源码需要替换的 Eigen 写法

核心源码中应优先消除以下 Eigen 专属写法：

### `iLQR/HessianCorrection.hpp`

将：

```cpp
matrix.diagonal().array() += minEigenvalue;
```

改为显式循环：

```cpp
for (int i = 0; i < Dimisions; ++i) {
  matrix(i, i) += minEigenvalue;
}
```

### `Dynamics/SystemDynamicsLinearizer.hpp`

将：

```cpp
Eigen::NumTraits<Scalar>::epsilon()
```

改为：

```cpp
std::numeric_limits<Scalar>::epsilon()
```

并引入：

```cpp
#include <limits>
```

### `Integration/RungeKuttaDormandPrince5.hpp`

将 `.array()` 与 `lpNorm<Eigen::Infinity>()` 替换为逐元素计算：

```cpp
Scalar maxError = Scalar(0);
for (int i = 0; i < XDim; ++i) {
  const Scalar scale =
      absTol + relTol * (std::abs(x_old(i)) +
                         std::abs(dt) * std::abs(dxdt_old(i)));
  maxError = std::max(maxError, std::abs(x_err(i) / scale));
}
return maxError;
```

### `Controller/LinearController.hpp`

将：

```cpp
uff.noalias() += k * x;
```

改为：

```cpp
uff += k * x;
```

### `Models/ThrustVector.hpp`

如果自定义矩阵补齐了 `head`、`tail`、`topRows`、`bottomRows`、`topLeftCorner`，这里可以保持大部分现有写法。

如果不补 Eigen-like API，则统一改成已有 `slice<P, Q>()` 风格。

## CMake 建议

主库目标不再公共链接 Eigen：

```cmake
target_link_libraries(iLQR_Include INTERFACE)
```

测试目标继续链接 Eigen，并启用兼容层：

```cmake
target_link_libraries(${test_name}
    PRIVATE
    GTest::gtest_main
    Eigen3::Eigen
    iLQR_Include
)

target_compile_definitions(${test_name}
    PRIVATE
    ILQR_ENABLE_EIGEN_COMPAT
)
```

如果迁移过程中需要快速回退测试到 Eigen 后端，可临时加入：

```cmake
target_compile_definitions(${test_name}
    PRIVATE
    ILQR_LINALG_BACKEND_EIGEN
)
```

但该方式不会覆盖自定义矩阵后端，因此只适合作为短期回退手段。

## 推荐迁移顺序

1. 在 `Types.hpp` 中加入后端选择机制，但先保留 Eigen 后端作为回退。
2. 为自定义 `Matrix` 和 `Vector` 补充最小 Eigen-like API。
3. 增加 `ILQR_ENABLE_EIGEN_COMPAT` 下的 Eigen 固定尺寸互转能力。
4. 修改 `Tests` 之外的 Eigen 专属写法。
5. 主库目标移除公共 Eigen 链接，测试目标继续链接 Eigen。
6. 编译非测试目标，确认核心源码不直接依赖 Eigen。
7. 编译并运行测试，确认 Eigen 测试代码可以通过兼容层与核心源码交互。

## 验证建议

建议分阶段验证：

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug --target iLQR
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

核心通过标准：

- `Tests` 之外不再直接包含 Eigen 头文件。
- 主库目标不再公共链接 `Eigen3::Eigen`。
- 固定尺寸核心求解路径使用自定义矩阵。
- `Tests` 中动态 QP 工具仍可使用 Eigen。
- 端到端测试数值行为与迁移前保持一致。

## 风险与约束

- 自定义矩阵目前是固定尺寸矩阵库，不应为了测试工具引入动态尺寸能力。
- `Eigen::Dynamic`、SVD、LU、动态 KKT 求解应暂时留在测试辅助代码中。
- Eigen-like API 只补项目真实需要的最小集合，不建议复刻完整 Eigen 接口。
- 迁移过程中应避免同时重构算法逻辑，确保数值差异可以归因到矩阵后端替换。
