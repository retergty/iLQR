# 离散状态空间适配方案

## 背景

当前工程的 iLQR 核心求解过程是离散时间形式。Riccati 反向递推使用的是离散线性二次模型：

$$
\delta x_{k+1} = A_k \delta x_k + B_k \delta u_k
$$

但当前模型层接口是连续时间动力学：

$$
\dot{x} = f(t, x, u)
$$

因此工程中存在一层连续到离散的转换：

1. 用户实现 `SystemDynamicsBase::computeFlowMap()`，返回状态导数。
2. `TimeTriggeredRollout` 通过 RK45 积分生成名义轨迹。
3. `SensitivityIntegrator` 将连续时间线性化离散化为 `ModelData`。
4. `DiscreteTimeRiccatiEquations` 对离散模型执行 backward pass。

对于本身就是采样系统的问题，例如力增量、上一拍控制量记忆、数字控制器状态更新等，直接写成离散状态空间更自然：

$$
x_{k+1} = f_d(t_k, x_k, u_k, \Delta t)
$$

目标是在不重写现有测试、不破坏连续时间模型的前提下，新增离散状态空间支持。

## 总体原则

不要把离散状态转移塞进现有 `computeFlowMap()` 接口。该接口语义是 `dx/dt`，如果返回 `x_{k+1}`，当前 rollout 和灵敏度离散化会把它再次积分，导致动力学错误。

推荐采用连续和离散两条并行路径：

```text
连续模型:
SystemDynamicsBase -> TimeTriggeredRollout -> SensitivityIntegrator -> ModelData -> Riccati

离散模型:
DiscreteSystemDynamicsBase -> DiscreteTimeRollout -> direct linearApproximation -> ModelData -> Riccati
```

Riccati 核心不需要改，因为它已经是离散时间递推。

## 建议增加转录层

为了避免 `iLQR.hpp` 同时承担连续动力学积分、连续线性化离散化和核心优化职责，建议增加一个代码层级，将“连续/离散模型如何变成离散 LQ 子问题”从底层 iLQR 优化流程中分离出来。

推荐分层：

```text
用户动力学层:
ContinuousDynamics / DiscreteDynamics

转录层:
ContinuousTranscription: RK rollout + sensitivityDiscretize + cost * dt
DiscreteTranscription: direct rollout + direct A/B + stage cost

核心优化层:
iLQRCore: 只处理离散 ModelData、离散 Riccati、控制器更新和线搜索
```

也就是说，连续系统和离散系统的差异应在转录层结束。进入 Riccati 后，核心求解器只看到统一的离散形式：

$$
\delta x_{k+1} = A_k \delta x_k + B_k \delta u_k
$$

### 推荐目录结构

```text
Dynamics/
  SystemDynamicsBase.hpp              // 现有连续 dx/dt 接口
  DiscreteSystemDynamicsBase.hpp      // 新增离散 x_next 接口

Transcription/
  ContinuousTranscription.hpp         // 连续模型转录为离散 LQ
  DiscreteTranscription.hpp           // 离散模型直接生成离散 LQ
  TranscriptionTraits.hpp             // 按 DynamicsMode 编译期选择

Rollout/
  TimeTriggeredRollout.hpp            // 现有连续 rollout
  DiscreteTimeRollout.hpp             // 新增离散 rollout
```

`iLQR` 内部不应直接关心动力学到底是连续还是离散，而是依赖转录器：

```cpp
using Transcription_t = typename Types::Transcription_t;

Transcription_t transcription_;
```

核心流程中只调用统一接口：

```cpp
transcription_.rollout(initTime, initState, finalTime, controller,
                       trajectory);

transcription_.approximateIntermediateLQ(
    optimalControlProblem_,
    targetTrajectory_,
    timeTrajectory[k],
    stateTrajectory[k],
    inputTrajectory[k],
    multiplierTrajectory[k],
    dt,
    modelDataTrajectory[k]);
```

### 连续转录器职责

连续转录器负责把连续系统转换成当前 Riccati 可用的离散模型：

```cpp
template <typename Descriptor>
class ContinuousTranscription {
 public:
  using Types = iLQRTypes<Descriptor>;

  int rollout(...);

  void approximateIntermediateLQ(..., Scalar dt, ModelData_t& modelData) {
    ModelData_t continuousModelData =
        Approximator_t::approximateIntermediateLQ(...);

    modelData.dynamics = discretizer_.sensitivityDiscretize(
        *problem.dynamicsPtr, time, state, input, dt);

    modelData.dynamics.f.setZero();
    modelData.cost = continuousModelData.cost * dt;
  }

 private:
  TimeTriggeredRollout_t rollout_;
  EK2DynamicsDiscretizer_t discretizer_;
};
```

这里的 `cost * dt` 只属于连续 running cost 的转录逻辑，不应散落在 iLQR 核心中。

### 离散转录器职责

离散转录器不做积分，也不做连续灵敏度离散化：

```cpp
template <typename Descriptor>
class DiscreteTranscription {
 public:
  using Types = iLQRTypes<Descriptor>;

  int rollout(...);

  void approximateIntermediateLQ(..., Scalar dt, ModelData_t& modelData) {
    ModelData_t stageModelData =
        Approximator_t::approximateIntermediateLQ(...);

    modelData.time = stageModelData.time;
    modelData.dynamics =
        problem.dynamicsPtr->linearApproximation(time, state, input, dt);

    modelData.dynamics.f.setZero();
    modelData.cost = stageModelData.cost;
  }

 private:
  DiscreteTimeRollout_t rollout_;
};
```

离散模式下，`modelData.cost` 默认表示每个采样步的 stage cost，因此不再乘 `dt`。

### 编译期选择

项目面向嵌入式和静态尺寸设计，推荐用模板 traits 在编译期选择转录器，不建议用运行时虚基类做主路径分发：

```cpp
struct ContinuousDynamics {};
struct DiscreteDynamics {};

template <typename Descriptor>
struct TranscriptionSelector;

template <typename Descriptor>
struct TranscriptionSelectorContinuous {
  using type = ContinuousTranscription<Descriptor>;
};

template <typename Descriptor>
struct TranscriptionSelectorDiscrete {
  using type = DiscreteTranscription<Descriptor>;
};
```

这样可以保持现有连续模型默认走 `ContinuousTranscription`，新离散模型显式走 `DiscreteTranscription`，并避免在核心求解路径中引入额外虚调用。

### 架构收益

增加转录层后，`iLQR` 的职责会收敛为：

```text
给定离散 ModelData
执行 Riccati backward pass
计算控制器修正
执行候选 rollout
比较 merit 并更新解
```

连续模型的 RK 积分、连续线性化离散化、`cost * dt` 缩放都属于 `ContinuousTranscription`；离散模型的直接状态递推、直接离散 A/B、离散 stage cost 都属于 `DiscreteTranscription`。这种分离更符合算法边界，也方便后续做嵌入式裁剪和性能优化。

## 修改范围

### 1. 增加动力学模式标记

在 `iLQRDescriptor.hpp` 中增加动力学模式 tag：

```cpp
struct ContinuousDynamics {};
struct DiscreteDynamics {};
```

扩展 `TranscriptionConfig`，增加第三个模板参数并默认使用连续模式：

```cpp
template <typename Dims_, typename Horizon_,
          typename DynamicsMode_ = ContinuousDynamics>
struct TranscriptionConfig {
  using Dims = Dims_;
  using Horizon = Horizon_;
  using DynamicsMode = DynamicsMode_;

  static constexpr int XDim = Dims::XDim;
  static constexpr int UDim = Dims::UDim;
  static constexpr std::size_t PredictLength = Horizon::PredictLength;
};
```

这样现有代码中的：

```cpp
using Config = TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<N>>;
```

仍然默认走连续路径，现有 Tests 不需要重写。

离散模型使用：

```cpp
using Config =
    TranscriptionConfig<Dimensions<XDim, UDim>, Horizon<N>, DiscreteDynamics>;
```

### 2. 增加离散动力学基类

新增 `Dynamics/DiscreteSystemDynamicsBase.hpp`。

接口语义应明确表示离散状态转移：

```cpp
template <typename Scalar, int XDim, int UDim>
class DiscreteSystemDynamicsBase {
 public:
  virtual ~DiscreteSystemDynamicsBase() = default;

  virtual Vector<Scalar, XDim> computeMap(
      Scalar t,
      const Vector<Scalar, XDim>& x,
      const Vector<Scalar, UDim>& u,
      Scalar dt) const = 0;

  virtual VectorFunctionLinearApproximation<Scalar, XDim, XDim, UDim>
  linearApproximation(
      Scalar t,
      const Vector<Scalar, XDim>& x,
      const Vector<Scalar, UDim>& u,
      Scalar dt) = 0;
};
```

其中 `linearApproximation()` 返回的是离散系统的一阶近似：

$$
x_{k+1} \approx f_d(x_k, u_k) + A_k \delta x_k + B_k \delta u_k
$$

这里的 `dfdx`、`dfdu` 分别是：

$$
A_k = \frac{\partial f_d}{\partial x}, \quad
B_k = \frac{\partial f_d}{\partial u}
$$

它们不是连续系统中的 `df/dx`、`df/du`。

### 3. 让问题对象按模式选择动力学指针

当前 `OptimalControlProblem` 中固定使用：

```cpp
SystemDynamicsBase<Scalar, XDim, UDim>* dynamicsPtr;
```

需要通过 traits 根据 `DynamicsMode` 选择类型：

```cpp
template <typename Scalar, int XDim, int UDim, typename DynamicsMode>
struct DynamicsTraits;

template <typename Scalar, int XDim, int UDim>
struct DynamicsTraits<Scalar, XDim, UDim, ContinuousDynamics> {
  using DynamicsBase = SystemDynamicsBase<Scalar, XDim, UDim>;
};

template <typename Scalar, int XDim, int UDim>
struct DynamicsTraits<Scalar, XDim, UDim, DiscreteDynamics> {
  using DynamicsBase = DiscreteSystemDynamicsBase<Scalar, XDim, UDim>;
};
```

然后：

```cpp
using DynamicsBase_t =
    typename DynamicsTraits<Scalar, XDim, UDim, DynamicsMode>::DynamicsBase;

DynamicsBase_t* dynamicsPtr = nullptr;
```

连续模式下类型与当前保持一致，已有模型和测试不会感知变化。

### 4. 增加离散 rollout

新增 `Rollout/DiscreteTimeRollout.hpp`，继承现有 `RolloutBase`，保持外部接口不变。

核心逻辑：

```cpp
for (size_t i = 0; i < numSteps + 1; ++i) {
  const Vector<Scalar, UDim> input =
      controller->computeInput(t, state);

  trajectory.timeTrajectory[i] = t;
  trajectory.stateTrajectory[i] = state;
  trajectory.inputTrajectory[i] = input;

  state = dynamics.computeMap(t, state, input, timeStep);
  t += timeStep;
}
```

注意最后一个点的输入可以继续按当前工程习惯重建，用于保持 `PredictLength + 1` 个输入槽位有值。

### 5. 按模式选择 rollout 类型

在 `iLQRTypes.hpp` 中增加 rollout traits：

```cpp
template <typename Scalar, int XDim, int UDim, typename DynamicsMode>
struct RolloutSelector;

template <typename Scalar, int XDim, int UDim>
struct RolloutSelector<Scalar, XDim, UDim, ContinuousDynamics> {
  using type = TimeTriggeredRollout<Scalar, XDim, UDim>;
};

template <typename Scalar, int XDim, int UDim>
struct RolloutSelector<Scalar, XDim, UDim, DiscreteDynamics> {
  using type = DiscreteTimeRollout<Scalar, XDim, UDim>;
};
```

然后令：

```cpp
using DynamicsBase_t = typename Traits::DynamicsBase_t;
using Rollout_t =
    typename RolloutSelector<Scalar, XDim, UDim, DynamicsMode>::type;
```

`iLQR` 内部成员从固定 `TimeTriggeredRollout_t rollout_;` 改成模式选择后的 `Rollout_t rollout_;`。

### 6. 修改 LQ 近似路径

当前 `iLQR::approximateIntermediateLQ()` 先生成连续 LQ 近似，再调用：

```cpp
modelData.dynamics = discretizer_.sensitivityDiscretize(
    system, time, state, input, timeStep);

modelData.cost = continuousTimeModelData.cost * timeStep;
```

离散模式不应进入 `SensitivityIntegrator`，也不应默认将代价再乘一次 `dt`。

建议拆成两个 worker：

```cpp
void lqWorkerContinuous(...);
void lqWorkerDiscrete(...);
```

连续模式：

```cpp
modelData.dynamics = discretizer_.sensitivityDiscretize(
    system, time, state, input, timeStep);
modelData.dynamics.f.setZero();
modelData.cost = continuousTimeModelData.cost * timeStep;
```

离散模式：

```cpp
modelData.time = continuousTimeModelData.time;
modelData.dynamics =
    system.linearApproximation(time, state, input, timeStep);
modelData.dynamics.f.setZero();
modelData.cost = continuousTimeModelData.cost;
```

这里保留 `modelData.dynamics.f.setZero()` 是为了与当前 deviation dynamics 逻辑一致。后续如果要支持显式 affine defect，再统一梳理 `ModelData::dynamics.f` 的语义。

### 7. 代价缩放策略

连续模式保持现状：

$$
L_k \approx l(x_k,u_k) \Delta t
$$

离散模式默认认为用户提供的是每步 stage cost：

$$
L_k = l_k(x_k,u_k)
$$

因此不乘 `dt`。

如果未来希望离散模型仍使用连续时间积分意义的 running cost，可以在 descriptor 中进一步增加 cost scaling policy，例如：

```cpp
struct ContinuousCostScaling {};
struct DiscreteCostScaling {};
```

但第一阶段不建议引入，避免扩大改动范围。

## 示例：矢量力增广离散模型

状态：

$$
X_k =
\begin{bmatrix}
v_k^W \\
F_{k-1}^B
\end{bmatrix}
\in \mathbb{R}^6
$$

输入：

$$
U_k = F_k^B \in \mathbb{R}^3
$$

离散动力学：

$$
X_{k+1} =
\begin{bmatrix}
v_k^W + \Delta t \left(\frac{1}{m}R_B^W U_k - g e_3 \right) \\
U_k
\end{bmatrix}
$$

雅可比：

$$
f_x =
\begin{bmatrix}
I & 0 \\
0 & 0
\end{bmatrix},
\quad
f_u =
\begin{bmatrix}
\frac{\Delta t}{m}R_B^W \\
I
\end{bmatrix}
$$

该模型不需要 RK 积分，也不需要连续灵敏度离散化。

## 测试策略

现有 Tests 不重写。因为默认模式是 `ContinuousDynamics`，所有旧 descriptor 默认继续走连续路径。

新增测试建议：

1. `DiscreteTimeRolloutTest`
   - 构造简单离散线性系统。
   - 验证 `x_{k+1} = A x_k + B u_k` 被逐步递推。

2. `DiscreteDynamicsApproximationTest`
   - 验证离散 `linearApproximation()` 返回的 `A`、`B` 不被再次乘 `dt`。

3. `DiscreteILQRSmokeTest`
   - 使用低维离散线性系统和二次代价。
   - 验证求解器可以完成一次或多次迭代，并生成有限代价。

回归测试：

```bash
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

## 风险点

### 语义混用

最大风险是连续和离散接口复用同一个函数名或同一条处理路径，导致离散系统被再次积分或再次离散化。应通过独立基类和 traits 分支避免。

### 代价是否乘 dt

连续 running cost 应乘 `dt`，离散 stage cost 默认不乘。两者必须在代码中明确分支。

### `ModelData::dynamics.f` 的含义

当前主求解器把离散动力学仿射项清零，使用 deviation dynamics：

$$
\delta x_{k+1} = A_k \delta x_k + B_k \delta u_k
$$

离散路径第一阶段应保持这个行为，避免影响 Riccati 和线搜索逻辑。

### 初始化器路径

`InitializerRollout` 本身已经是按步调用 `Initializer::compute()` 生成 `nextState`，不依赖连续积分，因此不需要为离散模型重写。

## 推荐实施顺序

1. 增加 `ContinuousDynamics` / `DiscreteDynamics` tag，并保证旧 descriptor 编译不变。
2. 增加 `DiscreteSystemDynamicsBase`。
3. 增加 `DiscreteTimeRollout`。
4. 在 `iLQRTypes` 中通过 traits 选择动力学基类和 rollout 类型。
5. 在 `iLQR` 的 LQ 近似阶段增加连续/离散分支。
6. 添加离散路径最小测试。
7. 运行现有测试，确认连续路径无回归。

## 结论

适配离散状态空间的工程难度中等，核心求解器不需要重写。最稳妥的方式是保留现有连续路径作为默认行为，通过 descriptor 显式开启离散路径。这样既能支持天然离散的模型，也能保持当前连续系统模型和测试集基本不变。
