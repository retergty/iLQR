# Augmented Lagrangian Vector Constraint Design

## Background

当前项目的增广拉格朗日约束实现是标量约束模型：

- `StateConstraint::getValue()` 返回 `Scalar`。
- `StateInputConstraint::getValue()` 返回 `Scalar`。
- `Multiplier` 保存一个标量 `lagrangian`。
- `LagrangianMetrics` 保存一个标量 `constraint`。
- 一个 `StateAugmentedLagrangian` 或 `StateInputAugmentedLagrangian` 对应一个标量约束。
- 多个标量约束通过 `StateAugmentedLagrangianCollection` 或 `StateInputAugmentedLagrangianCollection` 管理。

这种设计可以通过“按行拆分”的方式表达 `Ax = b`，但不原生支持一个拉格朗日项内部包含向量约束。

长期目标是修改核心接口，使其同时支持：

1. 一个增广拉格朗日项内部包含多个约束，即 `h(x)` 或 `h(x,u)` 是向量。
2. 一个 `OptimalControlProblem` 中同一类约束可以包含多个增广拉格朗日项。

目标结构如下：

```text
OptimalControlProblem
  stateInputEqualityLagrangianCollection
    term 0: h0(x, u) in R^2
    term 1: h1(x, u) in R^3
    term 2: h2(x, u) in R^1
```

## Design Principle

将约束分为两层：

- **Constraint term**：一个增广拉格朗日项，内部可以是向量约束 `h_j in R^Cj`。
- **Constraint collection**：一个最优控制问题中某一类约束的多个增广拉格朗日项。

因此：

- 一个 `AugmentedLagrangian` 不再表示一个标量约束，而是表示一个向量约束组。
- 一个 `Collection` 表示多个增广拉格朗日项的集合。
- `ConstraintConfig` 不再只保存标量约束数量，而应该保存每类约束的 term 列表以及每个 term 的维度。

## ConstraintConfig

建议引入 `ConstraintTerm` 和 `ConstraintGroupLayout`：

```cpp
template <int CDim_>
struct ConstraintTerm {
  static constexpr int CDim = CDim_;
};

template <typename... Terms_>
struct ConstraintGroupLayout {
  using Terms = std::tuple<Terms_...>;

  static constexpr int NumTerms = sizeof...(Terms_);
  static constexpr int TotalDim = (Terms_::CDim + ... + 0);

  template <std::size_t I>
  using Term = std::tuple_element_t<I, Terms>;
};
```

然后 `ConstraintConfig` 表示六类约束：

```cpp
template <
    typename StateEq_,
    typename StateIneq_,
    typename StateInputEq_,
    typename StateInputIneq_,
    typename FinalStateEq_,
    typename FinalStateIneq_>
struct ConstraintConfig {
  using StateEq = StateEq_;
  using StateIneq = StateIneq_;
  using StateInputEq = StateInputEq_;
  using StateInputIneq = StateInputIneq_;
  using FinalStateEq = FinalStateEq_;
  using FinalStateIneq = FinalStateIneq_;
};
```

示例：有两个 state-input equality 拉格朗日项，维度分别是 2 和 1。

```cpp
using ConstraintConfig_t = ConstraintConfig<
    ConstraintGroupLayout<>,                         // state equality
    ConstraintGroupLayout<>,                         // state inequality
    ConstraintGroupLayout<ConstraintTerm<2>,
                          ConstraintTerm<1>>,        // state-input equality
    ConstraintGroupLayout<>,                         // state-input inequality
    ConstraintGroupLayout<>,                         // final state equality
    ConstraintGroupLayout<>                          // final state inequality
>;
```

这里的语义是：

- state-input equality 有两个增广拉格朗日项。
- 第一个项内部约束维度为 2。
- 第二个项内部约束维度为 1。
- 总约束维度为 3，但它不是 3 个独立的拉格朗日对象。

## Multiplier

当前 `Multiplier` 是标量版本：

```cpp
template <typename Scalar>
struct Multiplier {
  Scalar penalty;
  Scalar lagrangian;
};
```

建议改为向量版本：

```cpp
template <typename Scalar, int CDim>
struct Multiplier {
  Scalar penalty = 1.0;
  Vector<Scalar, CDim> lagrangian = Vector<Scalar, CDim>::Zero();

  void setZero() {
    penalty = 1.0;
    lagrangian.setZero();
  }
};
```

如果需要支持每个约束分量不同的 penalty，也可以进一步改为：

```cpp
template <typename Scalar, int CDim>
struct Multiplier {
  Vector<Scalar, CDim> penalty = Vector<Scalar, CDim>::Ones();
  Vector<Scalar, CDim> lagrangian = Vector<Scalar, CDim>::Zero();
};
```

但第一阶段建议保留一个标量 `penalty`，保持数学和实现更简单。

## LagrangianMetrics

当前 `LagrangianMetrics` 保存标量约束值。建议改为：

```cpp
template <typename Scalar, int CDim>
struct LagrangianMetrics {
  Scalar penalty = 0.0;
  Vector<Scalar, CDim> constraint = Vector<Scalar, CDim>::Zero();

  void setZero() {
    penalty = 0.0;
    constraint.setZero();
  }
};
```

一个 term 的 penalty 仍然是标量，因为它最终贡献到 cost；constraint 是向量，因为它表示该 term 的约束残差。

## Constraint Interfaces

### State Constraint

建议将 `StateConstraint` 改为带约束维度 `CDim`：

```cpp
template <typename Scalar, int XDim, int CDim>
class StateConstraint {
 public:
  explicit StateConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateConstraint() = default;

  constexpr ConstraintOrder getOrder() const { return order_; }

  virtual Vector<Scalar, CDim> getValue(
      Scalar time,
      const Vector<Scalar, XDim>& state) const = 0;

  virtual VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0>
  getLinearApproximation(
      Scalar time,
      const Vector<Scalar, XDim>& state) const = 0;

 private:
  ConstraintOrder order_;
};
```

### State-Input Constraint

```cpp
template <typename Scalar, int XDim, int UDim, int CDim>
class StateInputConstraint {
 public:
  explicit StateInputConstraint(ConstraintOrder order) : order_(order) {}
  virtual ~StateInputConstraint() = default;

  constexpr ConstraintOrder getOrder() const { return order_; }

  virtual Vector<Scalar, CDim> getValue(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input) const = 0;

  virtual VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>
  getLinearApproximation(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input) const = 0;

 private:
  ConstraintOrder order_;
};
```

项目中已有 `VectorFunctionLinearApproximation`，因此线性向量约束可以直接复用。

## Penalty Interface

`AugmentedPenaltyBase` 应从标量形式扩展为向量形式：

```cpp
template <typename Scalar, int CDim>
class AugmentedPenaltyBase {
 public:
  virtual ~AugmentedPenaltyBase() = default;

  virtual Scalar getValue(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const = 0;

  virtual Vector<Scalar, CDim> getDerivative(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const = 0;

  virtual Matrix<Scalar, CDim, CDim> getSecondDerivative(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const = 0;

  virtual Vector<Scalar, CDim> updateMultiplier(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const = 0;

  virtual Vector<Scalar, CDim> initializeMultiplier() const = 0;
};
```

### Vector Quadratic Penalty

对于等式约束 `h = 0`，可使用：

```cpp
p(h) = -lambda^T h + 0.5 * rho * h^T h
```

实现建议：

```cpp
template <typename Scalar, int CDim>
class QuadraticPenalty final : public AugmentedPenaltyBase<Scalar, CDim> {
 public:
  struct Config {
    Scalar scale = 100.0;
    Scalar stepSize = 0.0;
  };

  explicit QuadraticPenalty(const Config& config) : config_(config) {}

  Scalar getValue(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const override {
    (void)time;
    return -lambda.dot(h) + Scalar(0.5) * config_.scale * h.squaredNorm();
  }

  Vector<Scalar, CDim> getDerivative(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const override {
    (void)time;
    return -lambda + config_.scale * h;
  }

  Matrix<Scalar, CDim, CDim> getSecondDerivative(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const override {
    (void)time;
    (void)lambda;
    (void)h;
    return config_.scale * Matrix<Scalar, CDim, CDim>::Identity();
  }

  Vector<Scalar, CDim> updateMultiplier(
      Scalar time,
      const Vector<Scalar, CDim>& lambda,
      const Vector<Scalar, CDim>& h) const override {
    (void)time;
    return lambda - config_.stepSize * config_.scale * h;
  }

  Vector<Scalar, CDim> initializeMultiplier() const override {
    return Vector<Scalar, CDim>::Zero();
  }

 private:
  Config config_;
};
```

## Penalty Chain Rule

`Penalty.hpp` 需要从标量链式法则改为向量链式法则。

对于线性向量约束：

```cpp
h = h0 + Jx * dx + Ju * du
```

令：

```cpp
g = dp / dh       // CDim
H = d2p / dh2     // CDim x CDim
```

则增广拉格朗日对标量 cost 的二次近似为：

```cpp
approx.f = p
approx.dfdx = Jx.transpose() * g
approx.dfdxx = Jx.transpose() * H * Jx
```

state-input 情况再加：

```cpp
approx.dfdu = Ju.transpose() * g
approx.dfduu = Ju.transpose() * H * Ju
approx.dfdux = Ju.transpose() * H * Jx
```

这一步是整个重构中最关键的数学连接：约束是向量，但最终仍然贡献一个标量 cost approximation。

## AugmentedLagrangian Term

一个 state-input 增广拉格朗日项建议定义为：

```cpp
template <typename Scalar, int XDim, int UDim, int CDim>
class StateInputAugmentedLagrangianInterface {
 public:
  virtual ~StateInputAugmentedLagrangianInterface() = default;

  virtual LagrangianMetrics<Scalar, CDim> getValue(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& multiplier) const = 0;

  virtual ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
  getQuadraticApproximation(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& multiplier) const = 0;

  virtual std::pair<Multiplier<Scalar, CDim>, Scalar> updateLagrangian(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Vector<Scalar, CDim>& constraint,
      const Multiplier<Scalar, CDim>& multiplier) const = 0;

  virtual Multiplier<Scalar, CDim> initializeLagrangian(Scalar time) const = 0;
};
```

实现类：

```cpp
template <typename Scalar, int XDim, int UDim, int CDim>
class StateInputAugmentedLagrangian final
    : public StateInputAugmentedLagrangianInterface<
          Scalar, XDim, UDim, CDim> {
 public:
  StateInputAugmentedLagrangian(
      StateInputConstraint<Scalar, XDim, UDim, CDim>* constraint,
      AugmentedPenaltyBase<Scalar, CDim>* penalty)
      : constraint_(constraint), penalty_(penalty) {}

  LagrangianMetrics<Scalar, CDim> getValue(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input,
      const Multiplier<Scalar, CDim>& multiplier) const override {
    const auto h = constraint_->getValue(time, state, input);
    const Scalar p =
        multiplier.penalty *
        penalty_->getValue(time, multiplier.lagrangian, h);
    return {p, h};
  }

 private:
  StateInputConstraint<Scalar, XDim, UDim, CDim>* constraint_;
  AugmentedPenaltyBase<Scalar, CDim>* penalty_;
};
```

`StateAugmentedLagrangian` 同理，只是不包含输入。

## Collection Structure

因为不同 term 的 `CDim` 可能不同，`StateInputAugmentedLagrangianCollection` 不适合继续使用同质 `std::array`。

建议改为 `std::tuple`：

```cpp
template <typename Scalar, int XDim, int UDim, typename GroupLayout>
class StateInputAugmentedLagrangianCollection;

template <typename Scalar, int XDim, int UDim, typename... Terms>
class StateInputAugmentedLagrangianCollection<
    Scalar, XDim, UDim, ConstraintGroupLayout<Terms...>> {
 public:
  using Layout = ConstraintGroupLayout<Terms...>;

  template <std::size_t I>
  using TermLayout = typename Layout::template Term<I>;

  template <std::size_t I>
  using TermPtr =
      const StateInputAugmentedLagrangianInterface<
          Scalar, XDim, UDim, TermLayout<I>::CDim>*;

  template <std::size_t I>
  void set(TermPtr<I> term) {
    std::get<I>(terms_) = term;
  }

 private:
  std::tuple<
      const StateInputAugmentedLagrangianInterface<
          Scalar, XDim, UDim, Terms::CDim>*...>
      terms_;
};
```

这里建议使用 `set<I>()` 而不是 runtime `add()`。

原因：

- `add()` 在运行时不知道下一个槽位的 `CDim`。
- `set<I>()` 可以在编译期检查 term 类型和维度是否匹配。
- 这符合当前项目固定维度、模板化、Eigen fixed-size 的风格。

注册示例：

```cpp
problem.equalityLagrangian.template set<0>(&contactLagrangian);
problem.equalityLagrangian.template set<1>(&speedLagrangian);
```

## Collection Quadratic Approximation

collection 的二次近似是所有 term 的 scalar cost approximation 求和：

```cpp
template <typename MultipliersTuple, std::size_t... Is>
ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
getQuadraticApproximationImpl(
    Scalar time,
    const Vector<Scalar, XDim>& state,
    const Vector<Scalar, UDim>& input,
    const MultipliersTuple& multipliers,
    std::index_sequence<Is...>) const {
  ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim> out;
  out.setZero();

  ((out += std::get<Is>(terms_)->getQuadraticApproximation(
        time, state, input, std::get<Is>(multipliers))), ...);

  return out;
}
```

外部接口：

```cpp
template <typename MultipliersTuple>
ScalarFunctionQuadraticApproximation<Scalar, XDim, UDim>
getQuadraticApproximation(
    Scalar time,
    const Vector<Scalar, XDim>& state,
    const Vector<Scalar, UDim>& input,
    const MultipliersTuple& multipliers) const {
  return getQuadraticApproximationImpl(
      time, state, input, multipliers,
      std::make_index_sequence<sizeof...(Terms)>{});
}
```

`getValue()`、`initializeLagrangian()`、`updateLagrangian()` 使用同样的 tuple 遍历模式。

## MultiplierCollection

对于每类约束，multipliers 应该是一个 tuple：

```cpp
template <typename Scalar, typename GroupLayout>
struct MultiplierGroup;

template <typename Scalar, typename... Terms>
struct MultiplierGroup<Scalar, ConstraintGroupLayout<Terms...>> {
  std::tuple<Multiplier<Scalar, Terms::CDim>...> terms;
};
```

整个单时刻 multiplier collection：

```cpp
template <typename Scalar, typename ConstraintConfig>
struct MultiplierCollection {
  MultiplierGroup<Scalar, typename ConstraintConfig::StateEq> stateEq;
  MultiplierGroup<Scalar, typename ConstraintConfig::StateIneq> stateIneq;
  MultiplierGroup<Scalar, typename ConstraintConfig::StateInputEq> stateInputEq;
  MultiplierGroup<Scalar, typename ConstraintConfig::StateInputIneq> stateInputIneq;
};
```

终端时刻可以复用同一结构，只是 state-input group 为空。

## Metrics

`Metrics` 同样应按 group 存 tuple：

```cpp
template <typename Scalar, typename GroupLayout>
struct LagrangianMetricsGroup;

template <typename Scalar, typename... Terms>
struct LagrangianMetricsGroup<Scalar, ConstraintGroupLayout<Terms...>> {
  std::tuple<LagrangianMetrics<Scalar, Terms::CDim>...> terms;
};
```

单时刻 metrics：

```cpp
template <typename Scalar, typename Dims, typename ConstraintConfig>
struct Metrics {
  Scalar cost;
  Vector<Scalar, Dims::XDim> dynamicsViolation;

  LagrangianMetricsGroup<Scalar, typename ConstraintConfig::StateEq> stateEqLagrangian;
  LagrangianMetricsGroup<Scalar, typename ConstraintConfig::StateIneq> stateIneqLagrangian;
  LagrangianMetricsGroup<Scalar, typename ConstraintConfig::StateInputEq> stateInputEqLagrangian;
  LagrangianMetricsGroup<Scalar, typename ConstraintConfig::StateInputIneq> stateInputIneqLagrangian;
};
```

`sumPenalties()` 需要改成 tuple fold：

```cpp
template <typename Scalar, typename... Terms>
Scalar sumPenalties(
    const LagrangianMetricsGroup<Scalar, ConstraintGroupLayout<Terms...>>& group) {
  return std::apply(
      [](const auto&... metrics) {
        return (Scalar(0) + ... + metrics.penalty);
      },
      group.terms);
}
```

## Example: Ax = b

原生向量 state equality constraint：

```cpp
template <typename Scalar, int XDim, int CDim>
class LinearStateEqualityConstraint final
    : public StateConstraint<Scalar, XDim, CDim> {
 public:
  using A_t = Matrix<Scalar, CDim, XDim>;
  using b_t = Vector<Scalar, CDim>;

  LinearStateEqualityConstraint(const A_t& A, const b_t& b)
      : StateConstraint<Scalar, XDim, CDim>(ConstraintOrder::Linear),
        A_(A),
        b_(b) {}

  Vector<Scalar, CDim> getValue(
      Scalar time,
      const Vector<Scalar, XDim>& state) const override {
    (void)time;
    return A_ * state - b_;
  }

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0>
  getLinearApproximation(
      Scalar time,
      const Vector<Scalar, XDim>& state) const override {
    VectorFunctionLinearApproximation<Scalar, CDim, XDim, 0> out;
    out.f = getValue(time, state);
    out.dfdx = A_;
    return out;
  }

 private:
  A_t A_;
  b_t b_;
};
```

配置示例：

```cpp
using ConstraintConfig_t = ConstraintConfig<
    ConstraintGroupLayout<ConstraintTerm<CDim>>,  // state equality: one vector term
    ConstraintGroupLayout<>,
    ConstraintGroupLayout<>,
    ConstraintGroupLayout<>,
    ConstraintGroupLayout<>,
    ConstraintGroupLayout<>>;
```

注册：

```cpp
static LinearStateEqualityConstraint<Scalar, XDim, CDim> constraint(A, b);
static QuadraticPenalty<Scalar, CDim> penalty(config);
static StateAugmentedLagrangian<Scalar, XDim, CDim> lagrangian(
    &constraint, &penalty);

problem.stateEqualityLagrangian.template set<0>(&lagrangian);
```

## Example: Ax + Bu = c

state-input equality constraint：

```cpp
template <typename Scalar, int XDim, int UDim, int CDim>
class LinearStateInputEqualityConstraint final
    : public StateInputConstraint<Scalar, XDim, UDim, CDim> {
 public:
  using A_t = Matrix<Scalar, CDim, XDim>;
  using B_t = Matrix<Scalar, CDim, UDim>;
  using c_t = Vector<Scalar, CDim>;

  LinearStateInputEqualityConstraint(
      const A_t& A,
      const B_t& B,
      const c_t& c)
      : StateInputConstraint<Scalar, XDim, UDim, CDim>(ConstraintOrder::Linear),
        A_(A),
        B_(B),
        c_(c) {}

  Vector<Scalar, CDim> getValue(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input) const override {
    (void)time;
    return A_ * state + B_ * input - c_;
  }

  VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim>
  getLinearApproximation(
      Scalar time,
      const Vector<Scalar, XDim>& state,
      const Vector<Scalar, UDim>& input) const override {
    VectorFunctionLinearApproximation<Scalar, CDim, XDim, UDim> out;
    out.f = getValue(time, state, input);
    out.dfdx = A_;
    out.dfdu = B_;
    return out;
  }

 private:
  A_t A_;
  B_t B_;
  c_t c_;
};
```

## Suggested Refactor Order

建议按以下顺序重构：

1. 新增 `ConstraintTerm`、`ConstraintGroupLayout`，重写 `ConstraintConfig`。
2. 修改 `Multiplier` 和 `LagrangianMetrics` 为带 `CDim` 的向量版本。
3. 修改 `StateConstraint` 和 `StateInputConstraint` 为向量约束接口。
4. 修改 `AugmentedPenaltyBase` 和 `QuadraticPenalty` 为向量版本。
5. 修改 `Penalty.hpp` 的链式法则，将向量约束映射为标量 cost 二次近似。
6. 修改 `StateAugmentedLagrangian` 和 `StateInputAugmentedLagrangian`，加入 `CDim` 模板参数。
7. 将 `StateAugmentedLagrangianCollection` 和 `StateInputAugmentedLagrangianCollection` 从 `std::array` 改为 `std::tuple`，使用 `set<I>()` 注册 term。
8. 修改 `MultiplierCollection`、`Metrics`、`ProblemMetrics`、`DualSolution`，让它们保存 tuple group。
9. 修改 `LinearQuadraticApproximator` 中对各类 lagrangian 的调用。
10. 修改 `OptimalControlProblemHelperFunction` 中初始化和更新 multiplier 的逻辑。
11. 修改模型示例，例如 `CircularKinematics`，让旧的单标量约束变成 `CDim = 1` 的向量约束。
12. 修改测试，先覆盖 `CDim = 1` 等价旧行为，再覆盖 `CDim > 1` 的 `Ax = b` 和 `Ax + Bu = c`。

## Recommendation

推荐采用编译期固定维度方案：

```cpp
ConstraintGroupLayout<ConstraintTerm<C0>, ConstraintTerm<C1>, ...>
```

collection 使用 `std::tuple`，通过 `set<I>()` 注册不同维度的增广拉格朗日项。

这个方案保留当前项目的模板化和 Eigen fixed-size 风格，同时实现：

- 一个增广拉格朗日项内部支持多个约束。
- 一个 `OptimalControlProblem` 中支持多个增广拉格朗日项。
- `Ax = b`、`Ax + Bu = c` 可以作为原生向量约束表达。

