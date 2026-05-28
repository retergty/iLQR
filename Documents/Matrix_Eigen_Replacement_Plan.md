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

## Eigen 代码替换方法

本章用于指导将固定尺寸 Eigen 写法改成自定义矩阵库写法。替换时应优先保留原变量名和矩阵维度，并显式补全 `matrix::Matrix` 或 `matrix::Vector` 的模板参数。

### 命名空间前缀替换

将 Eigen 类型迁移到自定义矩阵库时，先把 `Eigen::` 前缀替换为 `matrix::` 前缀，再根据目标类型补全自定义库的类型名与维度。

例如：

```cpp
Eigen::Matrix<double, 2, 2> A;
Eigen::Matrix<double, 3, 1> b;
```

替换为：

```cpp
matrix::Matrix<double, 2, 2> A;
matrix::Vector<double, 3> b;
```

注意 `Eigen::Matrix<Scalar, Rows, 1>` 应优先替换成 `matrix::Vector<Scalar, Rows>`，不要保留三模板参数的向量写法。

### `<<` 初始化优先替换为 `initializer_list`

Eigen 的逗号初始化依赖 `operator<<` 重载，自定义矩阵库不使用该写法。迁移时应优先改成 `initializer_list` 构造或赋值，并保持原来的行优先顺序；如果目标位置无法使用 `initializer_list`，再展开为逐元素赋值。

例如：

```cpp
Eigen::Matrix2d matrix;
matrix << 1.0, 2.0, 3.0, 4.0;
```

优先替换为：

```cpp
matrix::Matrix<double, 2, 2> matrix{1.0, 2.0, 3.0, 4.0};
```

也可以使用二维列表表达行列结构：

```cpp
matrix::Matrix<double, 2, 2> matrix{{1.0, 2.0}, {3.0, 4.0}};
```

如果对象已经存在，可使用 `initializer_list` 赋值：

```cpp
matrix = {1.0, 2.0, 3.0, 4.0};
```

如果当前位置无法使用 `initializer_list`，再展开为逐元素赋值：

```cpp
matrix(0, 0) = 1.0;
matrix(0, 1) = 2.0;
matrix(1, 0) = 3.0;
matrix(1, 1) = 4.0;
```

向量同样优先使用 `initializer_list`：

```cpp
Eigen::Vector3d v;
v << 1.0, 2.0, 3.0;
```

优先替换为：

```cpp
matrix::Vector<double, 3> v{1.0, 2.0, 3.0};
```

无法使用 `initializer_list` 时再逐项赋值：

```cpp
v(0) = 1.0;
v(1) = 2.0;
v(2) = 3.0;
```

### 固定尺寸别名手动展开

Eigen 的便捷别名需要手动展开为自定义矩阵类型，避免迁移后继续依赖 Eigen 类型定义。

常见替换如下：

```cpp
Eigen::Matrix2d matrix;
Eigen::Vector2d x;
Eigen::Matrix<double, XDim, UDim> B;
```

替换为：

```cpp
matrix::Matrix<double, 2, 2> matrix;
matrix::Vector<double, 2> x;
matrix::Matrix<double, XDim, UDim> B;
```

如果原类型是 `Eigen::Matrix<Scalar, Rows, 1>`，应展开为 `matrix::Vector<Scalar, Rows>`；如果列数不是 `1`，则展开为 `matrix::Matrix<Scalar, Rows, Cols>`。

### `Constant()` 替换为 `setAll()`

Eigen 的静态工厂函数 `Constant()` 需要替换为声明对象后调用 `setAll()`。由于当前自定义矩阵库的 `setAll()` 是成员函数且返回 `void`，不应直接写成链式表达式。

例如：

```cpp
auto v = Eigen::Vector<double, 1>::Constant(10.0);
```

替换为：

```cpp
matrix::Vector<double, 1> v;
v.setAll(10.0);
```

矩阵同理：

```cpp
auto W = Eigen::Matrix<double, XDim, XDim>::Constant(weight);
```

替换为：

```cpp
matrix::Matrix<double, XDim, XDim> W;
W.setAll(weight);
```

如果原代码必须在表达式内直接产生一个值，可用局部 lambda 包装初始化过程：

```cpp
const auto v = [] {
  matrix::Vector<double, 1> tmp;
  tmp.setAll(10.0);
  return tmp;
}();
```

### `Identity()` 替换为 `setIdentity()`

Eigen 的静态工厂函数 `Identity()` 需要替换为声明对象后调用 `setIdentity()`，并补全矩阵模板参数。

例如：

```cpp
auto Q = Matrix<Scalar, XDim, XDim>::Identity();
```

替换为：

```cpp
matrix::Matrix<Scalar, XDim, XDim> Q;
Q.setIdentity();
```

如果原代码已经使用项目级 `Matrix` 别名，可根据迁移范围选择保留别名或显式改成 `matrix::Matrix`。在直接替换 Eigen 代码时，推荐显式写出：

```cpp
matrix::Matrix<Scalar, XDim, XDim> Q;
Q.setIdentity();
```

表达式上下文中同样可使用 lambda：

```cpp
const auto Q = [] {
  matrix::Matrix<Scalar, XDim, XDim> tmp;
  tmp.setIdentity();
  return tmp;
}();
```

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
