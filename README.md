# iLQR 最优控制项目描述

## 项目定位

本项目是一个 iLQR 最优控制代码库，用于非线性系统的轨迹优化与约束优化求解。项目早期来源于 [OCS2](https://github.com/leggedrobotics/ocs2) 的代码体系和最优控制问题组织方式，但并不是对原框架的直接裁剪，而是在核心数据结构、求解流程、约束表达和内存管理策略上进行了重新设计。

整体而言，本项目目标是在更受控的工程假设下提供一个结构清晰、便于移植和分析的 iLQR 求解器实现。

## 设计目标与约束

本项目面向资源受限平台部署，设计上重点遵循以下约束：

- 禁止动态内存分配。核心求解流程不得依赖堆分配、运行时扩容或动态尺寸容器。
- 采用静态维度建模。状态维度、输入维度、预测步长和约束维度均通过模板参数或编译期配置表达。
- 使用固定容量数据结构。轨迹、控制器、乘子、模型近似、Riccati 缓存和性能指标等数据均使用预先确定容量的容器保存。

在约束处理方面，项目不采用复杂的严格约束求解框架，而是将等式和不等式约束建模为增广拉格朗日项。约束违反通过惩罚函数和拉格朗日乘子进入优化目标，并在迭代过程中更新乘子和惩罚参数，从而逐步改善约束满足程度。

## 算法依据

项目中的算法设计主要参考 `Documents` 目录下的内部说明文档。`Iterative_Linear_Quadratic_Regulator.md` 给出了基础 iLQR 方法，包括离散时间最优控制问题、贝尔曼最优性原理、动力学线性化、代价函数二次近似、反馈增益与前馈项求解，以及 value function 的反向递推关系。这部分内容对应代码中的线性二次近似、离散时间 Riccati 递推和控制律更新。

`AL_iLQR.md` 是约束处理部分的主要参考。该文档将增广拉格朗日方法引入 iLQR：内层在固定对偶变量和惩罚参数的条件下求解带增广项的 iLQR 子问题，外层根据轨迹约束违反情况更新乘子和惩罚参数。项目实现以这一思想为基础，将约束处理统一映射为可二次近似的代价项。

因此，本项目的整体算法结构可以概括为：以 iLQR 作为内层轨迹优化器，以增广拉格朗日方法作为外层约束协调机制，并通过静态尺寸数据结构实现确定性计算流程。

## 求解流程

核心求解器围绕名义轨迹反复执行以下阶段：

1. 根据初始状态和当前控制器生成或更新名义状态、输入和时间轨迹。
2. 根据动力学模式将名义轨迹附近的问题转录为离散 LQ 近似，并对目标函数和增广拉格朗日项进行二次近似。
3. 使用离散时间 Riccati 方程从终端时刻向初始时刻反向递推 value function。
4. 根据 Riccati 递推结果计算反馈增益和前馈控制修正。
5. 通过线搜索选择可接受的控制更新步长，并重新 rollout 得到候选轨迹。
6. 根据代价、约束违反和 merit 指标判断收敛情况。
7. 对增广拉格朗日乘子和惩罚参数进行更新，为下一轮约束优化提供对偶信息。

该流程将非线性约束最优控制问题分解为一系列固定结构的局部线性二次子问题，便于进行周期性执行、运行时间估计和内存预算。

## 模块结构

项目整体按功能划分模块。各模块通过编译期维度描述和固定长度轨迹容器连接，使状态维度、输入维度、预测步长和约束规模在类型层面保持明确。

### 核心求解器与类型系统

`iLQRCore.hpp` 是库用户推荐包含的公共聚合头，集中导出求解器、描述符、动力学、代价、约束、惩罚和基础线性代数类型。`iLQR/iLQR.hpp` 是主要求解器实现，负责组织名义轨迹初始化、问题近似、Riccati 反向递推、控制器更新、线搜索、性能指标计算和乘子更新等步骤。

`iLQR/iLQRDescriptor.hpp`、`iLQR/iLQRDescriptorTraits.hpp` 和 `iLQR/iLQRTypes.hpp` 构成静态类型描述体系。它们将标量类型、系统维度、预测长度、动力学模式和约束布局转换为统一的轨迹、控制器、乘子、模型数据和求解缓存类型，是项目实现固定尺寸和禁止动态内存分配的重要基础。`TranscriptionConfig` 支持 `ContinuousDynamics` 与 `DiscreteDynamics` 两种动力学模式，默认使用连续模式；约束布局采用 `ConstraintGroupLayout<ConstraintTerm<N>...>` 描述：每个 `ConstraintTerm<N>` 表示一个增广拉格朗日项内部的 `N` 维向量约束，同一类约束可以包含多个 term。

`iLQR/LinearAlgebraTypes.hpp` 提供核心线性代数入口，将项目内使用的 `Matrix<Scalar, Rows, Cols>` 和 `Vector<Scalar, Rows>` 别名绑定到自定义固定尺寸矩阵库。核心代码不再通过公共头依赖 Eigen，测试侧需要 Eigen 对照能力时通过 `Tests/Include` 中的转换工具接入。

`iLQR/DDPSetting.hpp`、`iLQR/DDPData.hpp` 和 `iLQR/HessianCorrection.hpp` 提供求解参数、迭代数据和 Hessian 修正逻辑，用于改善 Riccati 递推和控制更新的数值稳定性。求解器构造时需要外部提供 `OptimalControlProblem` 与初始化器对象，相关对象生命周期应长于求解器。

### 最优控制问题建模

`OptimalControl` 目录定义最优控制问题本体及辅助函数。`OptimalControlProblem.hpp` 将系统动力学、运行代价、终端代价和增广拉格朗日约束组织为统一问题对象，其中 `dynamicsPtr` 的基类由 `TranscriptionConfig::DynamicsMode` 在编译期决定；`OptimalControlProblemHelperFunction.hpp` 负责对偶解初始化、乘子更新和约束指标同步。

`ModelData` 目录保存求解流程中的节点模型数据、指标数据和乘子基础类型。`ModelData.hpp` 表示转录后的单节点离散 LQ 数据，`Metrics.hpp` 表示单节点代价、约束和拉格朗日项评估结果，`Multiplier.hpp` 表示增广拉格朗日乘子。

`OptimalControlData` 目录保存最优控制求解过程中的轨迹解、参考轨迹和全局性能数据。`PrimalSolution.hpp` 与 `DualSolution.hpp` 分别表示原始轨迹解与对偶乘子解，`Reference.hpp` 表示参考轨迹，`PerformanceIndex.hpp`、`PerformanceIndexEvaluator.hpp`、`ProblemMetrics.hpp` 和 `LagrangianMetrics.hpp` 用于记录和汇总代价、约束违反、merit 指标和拉格朗日项评估结果，为收敛判断和线搜索提供量化依据。连续模式下 performance 按时间戳做梯形积分，离散模式下按 stage 指标直接求和，以保持与转录层中的代价缩放语义一致。

### 线性二次近似

`Approximation` 目录负责函数近似格式和目标函数局部二次化。`LinearApproximation.hpp` 和 `QuadraticApproximation.hpp` 定义一阶与二阶近似的数据格式，`LinearQuadraticApproximator.hpp` 负责在给定节点上计算运行代价、状态代价和增广拉格朗日项的二次近似。完整的离散 LQ `ModelData` 由 `Transcription` 层组合目标函数近似与连续/离散动力学近似后生成。

该模块为转录层提供目标函数二次近似，其结果会与离散动力学线性化共同决定反馈增益、前馈修正和 value function 递推的计算基础。

### 转录层

`Transcription` 目录负责将不同动力学模式下的最优控制问题统一转录为 Riccati 递推所需的离散 LQ 节点。`ContinuousTranscription.hpp` 面向连续动力学：通过 `SensitivityIntegrator` 离散化连续线性化结果，并将 running cost 按时间步长缩放。`DiscreteTranscription.hpp` 面向离散动力学：直接调用离散系统的 `deviationLinearApproximation()` 生成偏差动力学，并将中间代价解释为离散 stage cost。`TranscriptionTraits.hpp` 根据 descriptor 中的 `DynamicsMode` 在编译期选择具体转录器。

### 动力学、Rollout 与积分

`Dynamics` 目录定义连续和离散动力学接口。`ControlledSystemBase.hpp` 与 `SystemDynamicsBase.hpp` 表达连续时间系统流映射及其线性化，`SystemDynamicsLinearizer.hpp` 和 `LinearSystemDynamics.hpp` 分别支持数值线性化和连续线性系统建模。`DiscreteSystemBase.hpp`、`DiscreteSystemDynamicsBase.hpp` 和 `DiscreteLinearSystemDynamics.hpp` 表达离散状态转移及其线性化。`DynamicsModeTraits.hpp` 根据 `ContinuousDynamics` 或 `DiscreteDynamics` 在编译期选择 `OptimalControlProblem::dynamicsPtr` 的基类类型。

`Rollout` 目录负责前向轨迹展开。`TimeTriggeredRollout.hpp` 根据当前控制器和固定时间步长积分连续系统状态，并可重建输入轨迹；`DiscreteTimeRollout.hpp` 直接调用离散状态转移逐节点推进状态；`RolloutTraits.hpp` 根据动力学模式在编译期选择 rollout 实现；`InitializerRollout.hpp` 用于在尚无有效控制器时生成初始名义轨迹。

`Integration` 目录提供连续系统积分工具，包括 Runge-Kutta Dormand-Prince 5 阶积分器、梯形积分、灵敏度积分和积分观察器等组件。

### 代价、约束与增广拉格朗日

`Cost` 目录定义运行代价、终端代价和代价集合，并向近似器提供代价值、一阶导数和二阶导数。

`Constraint` 目录定义状态约束和状态-输入约束接口，并通过约束阶次描述其可用的近似形式。约束以固定维度向量形式表达，例如 `StateInputConstraint<Scalar, XDim, UDim, CDim>` 返回 `Vector<Scalar, CDim>` 及对应的向量函数近似。约束不会作为独立硬约束直接进入 Riccati 求解，而是由增广拉格朗日模块转化为优化目标中的附加项。

`AugmentedLagrangian` 目录封装状态约束和状态-输入约束的增广拉格朗日表达。该模块将约束函数、惩罚函数和拉格朗日乘子组合起来，提供约束取值、二次近似和乘子更新接口，是本项目处理等式与不等式约束的核心机制。每个增广拉格朗日 term 可以绑定一个向量约束，乘子为同维度向量；多个 term 通过 collection 组织，并使用 `set<I>()` 在编译期按位置注册。

`Penalties` 目录提供增广拉格朗日所需的惩罚模型，包括二次惩罚、光滑绝对值惩罚、松弛铰链惩罚和 relaxed barrier 等形式。不同惩罚模型决定了约束违反在目标函数中的权重形式和乘子更新行为。`MultidimensionalPenalty` 负责将标量惩罚函数逐分量应用到向量约束，并通过链式法则生成最终标量代价的二次近似。

### Riccati 递推与搜索策略

`RiccatiEquations` 目录实现离散时间 Riccati 差分方程。`DiscreteTimeRiccatiEquations.hpp` 根据线性二次近似后的模型数据，从终端时刻反向递推 value function，并计算控制律中的反馈增益和前馈修正。`RiccatiModification.hpp` 描述递推过程中的正则化和修正项。

`SearchStrategy` 目录负责控制器更新后的步长选择。当前主要采用线搜索策略，在不同步长下执行 rollout，并根据 merit 指标和下降条件选择可接受的控制更新，以提高 iLQR 迭代稳定性。

### 控制器与初始化

`Controller` 目录定义控制器接口和线性反馈控制器。`LinearController.hpp` 表示时变线性控制律 `u(t, x) = u_ff(t) + K(t)x`，并支持基于时间插值的控制输入计算。

`Initialization` 目录定义轨迹初始化器接口及默认实现，用于在求解器尚未生成有效控制器时构造初始输入和状态轨迹。

### 模型、工具与测试

`ExampleModels` 目录包含用于验证和演示的任务模型，例如双积分器到达任务、圆形运动学模型和实验性模型。这些模型将动力学、代价和约束组合起来，用于端到端验证求解器行为。

`AutomaticDifferentation` 目录当前主要提供有限差分数值求导工具，用于在缺少解析导数时计算动力学或函数的 Jacobian。`Misc` 目录包含线性代数、数值工具和插值等通用辅助函数。`IntrusiveList` 目录提供侵入式链表结构，用于支持无需额外动态分配的容器组织方式。

`Tests/Include` 目录提供测试辅助头文件，包括复用主项目转录层的 QP 离散转录、QP 求解器以及与 OCS2 测试工具相关的对照实现。QP 工具提供连续动力学和离散动力学两组显式别名，例如 `QpContinuousDynamicsOptimalControlProblem_t` 与 `QpDiscreteDynamicsOptimalControlProblem_t`。`Tests` 目录覆盖近似器、动力学、积分器、Riccati 递推、惩罚函数、rollout、线性系统 iLQR、连续/离散 correctness 对照和端到端求解流程，用于验证各模块的数值正确性和集成行为。

`CMakeLists.txt` 负责组织项目构建、依赖配置和测试目标。`build` 与 `out` 等目录为本地构建产物，不属于核心源码模块。

## 使用方法

### 构建项目

项目使用 CMake 组织构建，要求 CMake 版本不低于 3.20，并使用 C++17 标准。核心库通过 `iLQR::iLQR` 接口目标暴露，使用方可包含根目录聚合头 `iLQRCore.hpp`。核心代码不依赖 Eigen3；测试和 QP 对照工具依赖 Eigen3 与 GoogleTest，测试构建流程会按需查找或获取这些依赖。

外部 CMake 工程接入时，推荐链接命名空间目标：

```cmake
target_link_libraries(your_target PRIVATE iLQR::iLQR)
```

推荐使用仓库中的 CMake Preset 进行构建。以 GCC Debug 配置为例：

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
```

Release 版本可使用：

```bash
cmake --preset gcc-release
cmake --build --preset gcc-release
```

构建完成后，生成文件位于 `out/build/<preset-name>` 目录中。测试目标由 `Tests/CMakeLists.txt` 自动生成，并通过 `ILQR_BUILD_TESTING` 控制是否参与构建：独立构建本项目时默认开启，作为 `add_subdirectory()` 子项目集成时默认关闭。开启测试后仍遵循 CMake 标准 `BUILD_TESTING` 选项。开发用编译和链接选项通过 `ILQR_ENABLE_DEV_OPTIONS` 控制，默认仅在独立构建时开启，并由测试目标显式链接 `iLQR::Options` 使用，避免影响父项目目标。性能敏感的核心模板代码可通过 `ILQR_FORCE_OPTIMIZATION` 向消费目标传播 `-O3` 和 `NDEBUG`：独立构建本项目时默认关闭以便运行调试测试，作为 `add_subdirectory()` 子项目集成时默认开启，防止上层工程误用 Debug 配置导致 iLQR 运行性能异常；如需逐步调试核心代码，可显式配置 `-DILQR_FORCE_OPTIMIZATION=OFF`。

### 运行测试

完成配置和构建后，可以通过 CTest 运行测试：

```bash
ctest --preset gcc-debug
```

也可以进入对应构建目录后直接执行单个测试程序，例如：

```bash
./out/build/gcc-debug/Tests/iLQREndToEndTest
```

测试用例覆盖了动力学、连续/离散 rollout、转录层、代价函数、Riccati 递推、线性二次近似、惩罚函数、QP 对照、连续/离散 correctness 和端到端 iLQR 求解过程。使用或修改核心模块后，应优先运行相关测试确认数值行为是否保持一致。

### 接入求解器

在新的最优控制任务中使用求解器时，通常需要完成以下步骤：

1. 定义问题维度和预测步长，构造 `iLQRDescriptor`。
2. 提供系统动力学对象。连续系统实现 `SystemDynamicsBase` 所需的流映射和线性化接口；离散系统实现 `DiscreteSystemDynamicsBase` 所需的状态转移和离散线性化接口。
3. 提供轨迹初始化器，用于生成初始名义轨迹。
4. 配置 `DDPSettings`，包括时间步长、最大迭代次数、收敛阈值和搜索策略等参数。
5. 构造 `OptimalControlProblem`，注册动力学、代价函数和约束项。
6. 构造 `iLQR<Descriptor>` 求解器。
7. 设置参考时间、状态和输入轨迹。
8. 调用 `run(initTime, initState)` 执行求解。
9. 通过 `primalSolution()` 读取优化后的状态轨迹、输入轨迹和线性反馈控制器。

典型的最小调用形式如下：

```cpp
#include "iLQRCore.hpp"

#include <array>
#include <cstddef>

using Scalar = double;
static constexpr int XDim = 2;
static constexpr int UDim = 2;
static constexpr std::size_t PredictLength = 5;

using Descriptor =
    iLQRDescriptor<Scalar,
                   TranscriptionConfig<Dimensions<XDim, UDim>,
                                       Horizon<PredictLength>>>;

using Solver = iLQR<Descriptor>;
using Problem = Solver::OptimalControlProblem_t;

const Matrix<Scalar, XDim, XDim> A = Matrix<Scalar, XDim, XDim>::Identity();
const Matrix<Scalar, XDim, UDim> B = Matrix<Scalar, XDim, UDim>::Identity();
LinearSystemDynamics<Scalar, XDim, UDim> dynamics(A, B);
DefaultInitializer<Scalar, XDim, UDim> initializer;
Problem problem;
problem.dynamicsPtr = &dynamics;

DDPSettings<Scalar> settings;
settings.timeStep_ = 0.01;
settings.maxNumIterations_ = 30;
settings.strategy_ = SearchStrategyType::LINE_SEARCH;

const Matrix<Scalar, XDim, XDim> Q = Matrix<Scalar, XDim, XDim>::Identity();
const Matrix<Scalar, UDim, UDim> R = Matrix<Scalar, UDim, UDim>::Identity();
const Matrix<Scalar, XDim, XDim> Qf =
    Scalar(10.0) * Matrix<Scalar, XDim, XDim>::Identity();
QuadraticStateInputCost<Scalar, XDim, UDim, PredictLength + 1> runningCost(
    Q, R, 0);
QuadraticStateCost<Scalar, XDim, PredictLength + 1> finalCost(Qf, 0);
problem.cost.add(runningCost);
problem.finalCost.add(finalCost);

Solver solver(settings, problem, &initializer);

std::array<Scalar, PredictLength + 1> timeTrajectory{};
std::array<Vector<Scalar, XDim>, PredictLength + 1> stateTrajectory{};
std::array<Vector<Scalar, UDim>, PredictLength + 1> inputTrajectory{};

// 填充参考轨迹后写入求解器。
solver.setDesireTrajectory(timeTrajectory, stateTrajectory, inputTrajectory);

Vector<Scalar, XDim> initState = Vector<Scalar, XDim>::Zero();
solver.run(0.0, initState);

const auto& solution = solver.primalSolution();
const auto& controller = solution.controller_;
```

在该示例中，`Dimensions<XDim, UDim>` 表示状态维度和输入维度，`Horizon<PredictLength>` 表示预测区间包含 `PredictLength` 个控制阶段和 `PredictLength + 1` 个时间节点。因此，与参考轨迹相关的数组长度应为 `PredictLength + 1`。由于 `TranscriptionConfig` 未显式指定第三个模板参数，该示例默认使用 `ContinuousDynamics`。

若模型本身已经是离散时间系统，可显式选择离散动力学模式：

```cpp
using DiscreteDescriptor =
    iLQRDescriptor<double,
                   TranscriptionConfig<Dimensions<2, 2>,
                                       Horizon<5>,
                                       DiscreteDynamics>>;

using DiscreteSolver = iLQR<DiscreteDescriptor>;
using DiscreteProblem = DiscreteSolver::OptimalControlProblem_t;

DiscreteLinearSystemDynamics<double, 2, 2> discreteDynamics(Ad, Bd);
DefaultInitializer<double, 2, 2> initializer;

DiscreteProblem problem;
problem.dynamicsPtr = &discreteDynamics;

DiscreteSolver solver(settings, problem, &initializer);
```

离散模式下，rollout 使用 `computeMap(t, x, u, dt)` 推进状态，中间代价默认按每步 stage cost 解释，不再自动乘以 `dt`。连续模式下，中间代价按 running cost density 解释，并由连续转录器乘以时间步长。

### 添加代价与约束

运行代价和终端代价通过 `OptimalControlProblem` 中的代价集合注册。对于二次型跟踪任务，可以使用 `QuadraticStateInputCost` 和 `QuadraticStateCost`；对于更复杂的任务，需要派生相应的代价接口，并提供代价值、一阶导数和二阶导数。

约束通过增广拉格朗日模块接入。状态约束和状态-输入约束应先实现对应的向量约束接口，再与惩罚函数组合为增广拉格朗日项，最后注册到 `OptimalControlProblem` 中的等式或不等式约束集合。求解器会在局部近似阶段将这些约束项转换为代价函数中的二次近似项，并在迭代过程中更新乘子。

约束配置采用两层结构：

```cpp
using ConstraintConfig_t = ConstraintConfig<
    StateConstraintConfig<ConstraintLayout<
        ConstraintGroupLayout<>,                         // state equality
        ConstraintGroupLayout<>>>,                       // state inequality
    StateInputConstraintConfig<ConstraintLayout<
        ConstraintGroupLayout<ConstraintTerm<1>>,         // state-input equality
        ConstraintGroupLayout<>>>,                       // state-input inequality
    FinalStateConstraintConfig<ConstraintLayout<
        ConstraintGroupLayout<>,                         // final state equality
        ConstraintGroupLayout<>>>>;                      // final state inequality
```

其中 `ConstraintGroupLayout<ConstraintTerm<1>>` 表示该类约束包含一个增广拉格朗日 term，term 内部约束维度为 1。如果需要一个 3 维约束项，可写成 `ConstraintGroupLayout<ConstraintTerm<3>>`；如果同一类约束有两个 term，维度分别为 2 和 1，则写成 `ConstraintGroupLayout<ConstraintTerm<2>, ConstraintTerm<1>>`。

以一个状态-输入等式约束为例：

```cpp
using Constraint_t = MyStateInputConstraint<double>;  // 派生自 StateInputConstraint<double, XDim, UDim, CDim>
using Penalty_t = QuadraticPenalty<double>;
using Lagrangian_t =
    StateInputAugmentedLagrangian<double, XDim, UDim, CDim>;

static Constraint_t constraint;
static Penalty_t penalty(typename Penalty_t::Config{100.0, 0.1});
static Lagrangian_t lagrangian(&constraint, &penalty);

problem.equalityLagrangian.template set<0>(&lagrangian);
```

注意，代价集合仍使用 `add()` 注册：

```cpp
problem.cost.add(runningCost);
problem.finalCost.add(finalCost);
```

而增广拉格朗日约束集合使用 `set<I>()` 注册，因为不同 term 的约束维度可能不同，需要在编译期确定位置和类型。

### 内存约束

使用或扩展本项目时，应保持问题规模在编译期确定。新增模块应优先使用固定尺寸 `Matrix` / `Vector`、`std::array` 和预分配缓存，避免在核心求解流程中引入动态内存分配。

### 注释与文档约定

源码注释统一采用中文描述，保留现有 Doxygen 标签、数学公式、文件名、类型名和代码标识符的原始写法。对于从 OCS2 继承或参考而来的模块，注释应说明当前实现中的职责和约束，不直接照搬外部框架语境。

新增或修改代码时，应优先保持原有注释风格：短行内注释用于解释局部实现意图，函数或类型级注释用于说明输入输出、数学含义和固定尺寸假设。涉及求解流程、约束布局、内存策略或模块职责变化时，应同步更新本文档。
