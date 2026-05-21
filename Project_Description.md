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
2. 在名义轨迹附近对系统动力学进行线性化，并对代价函数和增广拉格朗日项进行二次近似。
3. 使用离散时间 Riccati 方程从终端时刻向初始时刻反向递推 value function。
4. 根据 Riccati 递推结果计算反馈增益和前馈控制修正。
5. 通过线搜索选择可接受的控制更新步长，并重新 rollout 得到候选轨迹。
6. 根据代价、约束违反和 merit 指标判断收敛情况。
7. 对增广拉格朗日乘子和惩罚参数进行更新，为下一轮约束优化提供对偶信息。

该流程将非线性约束最优控制问题分解为一系列固定结构的局部线性二次子问题，便于进行周期性执行、运行时间估计和内存预算。

## 模块结构

项目整体按功能划分模块。各模块通过编译期维度描述和固定长度轨迹容器连接，使状态维度、输入维度、预测步长和约束规模在类型层面保持明确。

### 核心求解器与类型系统

`iLQR.hpp` 是主要求解器实现，负责组织名义轨迹初始化、问题近似、Riccati 反向递推、控制器更新、线搜索、性能指标计算和乘子更新等步骤。`iLQR.cpp` 提供最小化入口示例，用于展示求解器的基本构造方式。

`iLQRDescriptor.hpp`、`iLQRDescriptorTraits.hpp` 和 `iLQRTypes.hpp` 构成静态类型描述体系。它们将标量类型、系统维度、预测长度和约束布局转换为统一的轨迹、控制器、乘子、模型数据和求解缓存类型，是项目实现固定尺寸和禁止动态内存分配的重要基础。

`DDPSetting.hpp`、`DDPData.hpp` 和 `HessianCorrection.hpp` 提供求解参数、迭代数据和 Hessian 修正逻辑，用于改善 Riccati 递推和控制更新的数值稳定性。

### 最优控制问题建模

`OptimalControl` 目录定义最优控制问题及其求解状态。`OptimalControlProblem.hpp` 将系统动力学、运行代价、终端代价、参考轨迹和增广拉格朗日约束组织为统一问题对象。`PrimalSolution.hpp` 和 `DualSolution.hpp` 分别表示原始轨迹解与对偶乘子解。

`PerformanceIndex.hpp`、`Metrics.hpp`、`ProblemMetrics.hpp` 和 `LagrangianMetrics.hpp` 用于记录代价、约束违反、merit 指标和拉格朗日项评估结果，为收敛判断和线搜索提供量化依据。

### 线性二次近似

`Approximation` 目录负责将非线性最优控制问题转换为局部线性二次子问题。`LinearApproximation.hpp` 和 `QuadraticApproximation.hpp` 定义一阶与二阶近似的数据格式，`LinearQuadraticApproximator.hpp` 在每个时间节点上组合动力学线性化、代价二次化和增广拉格朗日项二次近似，生成 Riccati 递推所需的 `ModelData`。

该模块是非线性问题与 LQ 子问题之间的接口，其输出直接决定反馈增益、前馈修正和 value function 递推的计算基础。

### 动力学、Rollout 与积分

`Dynamics` 目录定义受控系统和动力学线性化接口。`ControlledSystemBase.hpp` 提供系统流映射，`SystemDynamicsBase.hpp` 增加线性化接口，`SystemDynamicsLinearizer.hpp` 和 `LinearSystemDynamics.hpp` 分别支持数值线性化和线性系统建模。

`Rollout` 目录负责前向轨迹展开。`TimeTriggeredRollout.hpp` 根据当前控制器和固定时间步长积分系统状态，并可重建输入轨迹；`InitializerRollout.hpp` 用于在尚无有效控制器时生成初始名义轨迹。

`Integration` 目录提供连续系统积分工具，包括 Runge-Kutta Dormand-Prince 5 阶积分器、梯形积分、灵敏度积分和积分观察器等组件。

### 代价、约束与增广拉格朗日

`Cost` 目录定义运行代价、终端代价和代价集合，并向近似器提供代价值、一阶导数和二阶导数。

`Constraint` 目录定义状态约束和状态-输入约束接口，并通过约束阶次描述其可用的近似形式。约束不会作为独立硬约束直接进入 Riccati 求解，而是由增广拉格朗日模块转化为优化目标中的附加项。

`AugmentedLagrangian` 目录封装状态约束和状态-输入约束的增广拉格朗日表达。该模块将约束函数、惩罚函数和拉格朗日乘子组合起来，提供约束取值、二次近似和乘子更新接口，是本项目处理等式与不等式约束的核心机制。

`Penalties` 目录提供增广拉格朗日所需的惩罚模型，包括二次惩罚、光滑绝对值惩罚、松弛铰链惩罚和 relaxed barrier 等形式。不同惩罚模型决定了约束违反在目标函数中的权重形式和乘子更新行为。

### Riccati 递推与搜索策略

`RiccatiEquations` 目录实现离散时间 Riccati 差分方程。`DiscreteTimeRiccatiEquations.hpp` 根据线性二次近似后的模型数据，从终端时刻反向递推 value function，并计算控制律中的反馈增益和前馈修正。`RiccatiModification.hpp` 描述递推过程中的正则化和修正项。

`SearchStrategy` 目录负责控制器更新后的步长选择。当前主要采用线搜索策略，在不同步长下执行 rollout，并根据 merit 指标和下降条件选择可接受的控制更新，以提高 iLQR 迭代稳定性。

### 控制器与初始化

`Controller` 目录定义控制器接口和线性反馈控制器。`LinearController.hpp` 表示时变线性控制律 `u(t, x) = u_ff(t) + K(t)x`，并支持基于时间插值的控制输入计算。

`Initialization` 目录定义轨迹初始化器接口及默认实现，用于在求解器尚未生成有效控制器时构造初始输入和状态轨迹。

### 模型、工具与测试

`Models` 目录包含用于验证和演示的任务模型，例如双积分器到达任务、圆形运动学模型和实验性模型。这些模型将动力学、代价和约束组合起来，用于端到端验证求解器行为。

`AutomaticDifferentation` 目录当前主要提供有限差分数值求导工具，用于在缺少解析导数时计算动力学或函数的 Jacobian。`Misc` 目录包含线性代数、数值工具和插值等通用辅助函数。`IntrusiveList` 目录提供侵入式链表结构，用于支持无需额外动态分配的容器组织方式。

`TestTools` 目录提供测试辅助工具，包括 QP 离散转录、QP 求解器以及与 OCS2 测试工具相关的对照实现。`Tests` 目录覆盖近似器、动力学、积分器、Riccati 递推、惩罚函数、rollout、线性系统 iLQR 和端到端求解流程，用于验证各模块的数值正确性和集成行为。

`CMakeLists.txt` 负责组织项目构建、依赖配置和测试目标。`build` 与 `out` 等目录为本地构建产物，不属于核心源码模块。

## 使用方法

### 构建项目

项目使用 CMake 组织构建，要求 CMake 版本不低于 3.20，并使用 C++17 标准。项目依赖 Eigen3；如果系统中未找到 Eigen3，构建脚本会通过 `FetchContent` 获取 Eigen。测试部分依赖 GoogleTest，同样由测试构建流程自动获取。

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

构建完成后，生成文件位于 `out/build/<preset-name>` 目录中。主示例目标为 `iLQR`，测试目标由 `Tests/CMakeLists.txt` 自动生成。

### 运行测试

完成配置和构建后，可以通过 CTest 运行测试：

```bash
ctest --preset gcc-debug
```

也可以进入对应构建目录后直接执行单个测试程序，例如：

```bash
./out/build/gcc-debug/Tests/iLQREndToEndTest
```

测试用例覆盖了动力学、积分器、rollout、代价函数、Riccati 递推、线性二次近似、惩罚函数和端到端 iLQR 求解过程。使用或修改核心模块后，应优先运行相关测试确认数值行为是否保持一致。

### 接入求解器

在新的最优控制任务中使用求解器时，通常需要完成以下步骤：

1. 定义问题维度和预测步长，构造 `iLQRDescriptor`。
2. 提供系统动力学对象，并实现 `SystemDynamicsBase` 所需的流映射和线性化接口。
3. 提供轨迹初始化器，用于生成初始名义轨迹。
4. 配置 `DDPSettings`，包括时间步长、最大迭代次数、收敛阈值和搜索策略等参数。
5. 构造 `iLQR<Descriptor>` 求解器，并向 `optimalControlProblem_` 注册代价函数和约束项。
6. 设置参考时间、状态和输入轨迹。
7. 调用 `run(initTime, initState)` 执行求解。
8. 通过 `primalSolution()` 读取优化后的状态轨迹、输入轨迹和线性反馈控制器。

典型的最小调用形式如下：

```cpp
using Descriptor =
    iLQRDescriptor<double,
                   TranscriptionConfig<Dimensions<2, 2>, Horizon<5>>>;

using Solver = iLQR<Descriptor>;

LinearSystemDynamics<double, 2, 2> dynamics(A, B);
DefaultInitializer<double, 2, 2> initializer;

DDPSettings<double> settings;
settings.timeStep_ = 0.01;
settings.maxNumIterations_ = 30;
settings.strategy_ = SearchStrategyType::LINE_SEARCH;

Solver solver(settings, &dynamics, &initializer);

QuadraticStateInputCost<double, 2, 2, 6> runningCost(Q, R);
QuadraticStateCost<double, 2, 6> finalCost(Qf);
solver.optimalControlProblem_.cost.add(runningCost);
solver.optimalControlProblem_.finalCost.add(finalCost);

std::array<double, 6> timeTrajectory;
std::array<Eigen::Vector2d, 6> stateTrajectory;
std::array<Eigen::Vector2d, 6> inputTrajectory;

// 填充参考轨迹后写入求解器。
solver.setDesireTrajectory(timeTrajectory, stateTrajectory, inputTrajectory);

Eigen::Vector2d initState;
solver.run(0.0, initState);

const auto& solution = solver.primalSolution();
const auto& controller = solution.controller_;
```

在该示例中，`Dimensions<2, 2>` 表示状态维度为 2、输入维度为 2，`Horizon<5>` 表示预测区间包含 5 个控制阶段和 6 个时间节点。因此，与轨迹相关的数组长度应为 `PredictLength + 1`。

### 添加代价与约束

运行代价和终端代价通过 `optimalControlProblem_` 中的代价集合注册。对于二次型跟踪任务，可以使用 `QuadraticStateInputCost` 和 `QuadraticStateCost`；对于更复杂的任务，需要派生相应的代价接口，并提供代价值、一阶导数和二阶导数。

约束通过增广拉格朗日模块接入。状态约束和状态-输入约束应先实现对应的约束接口，再与惩罚函数组合为增广拉格朗日项，最后注册到 `OptimalControlProblem` 中的等式或不等式约束集合。求解器会在局部近似阶段将这些约束项转换为代价函数中的二次近似项，并在迭代过程中更新乘子。

### 内存约束

使用或扩展本项目时，应保持问题规模在编译期确定。新增模块应优先使用固定尺寸 Eigen 类型、`std::array` 和预分配缓存，避免在核心求解流程中引入动态内存分配。
