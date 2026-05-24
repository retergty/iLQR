/**
 * @file RungeKuttaDormandPrince5.hpp
 * @brief Runge-Kutta Dormand-Prince 5(4) 单步积分器与固定步长积分封装。
 */
#pragma once

#include "Integration.hpp"

/**
 * @brief Dormand-Prince 5 阶单步器：给定 (t, x, dxdt) 与步长
 * dt，计算下一状态与导数。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 */
template <typename Scalar, int XDim>
class RungeKuttaDormandPrince5Stepper {
 public:
  RungeKuttaDormandPrince5Stepper() = default;
  /**
   * @brief 执行一步 Dormand-Prince 积分。
   * @param [in] system 系统动力学。
   * @param [in] x0 当前状态。
   * @param [in] dxdt 当前状态导数。
   * @param [in] t 当前时间。
   * @param [in] dt 步长。
   * @param [out] x_out 下一状态。
   * @param [out] dxdt_out 下一状态的导数。
   * @param [in] computeDxdtOut 是否计算下一状态导数；最后一步可跳过。
   */
  void doStep(const OdeBase<Scalar, XDim>& system,
              const Vector<Scalar, XDim>& x0, const Vector<Scalar, XDim>& dxdt,
              const Scalar t, const Scalar dt, Vector<Scalar, XDim>& x_out,
              Vector<Scalar, XDim>& dxdt_out, bool computeDxdtOut = true) {
    /* Runge-Kutta Dormand-Prince Butcher 表常数。
     * https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method */
    constexpr Scalar a2 = 1.0 / 5;
    constexpr Scalar a3 = 3.0 / 10;
    constexpr Scalar a4 = 4.0 / 5;
    constexpr Scalar a5 = 8.0 / 9;

    constexpr Scalar b21 = 1.0 / 5;

    constexpr Scalar b31 = 3.0 / 40;
    constexpr Scalar b32 = 9.0 / 40;

    constexpr Scalar b41 = 44.0 / 45;
    constexpr Scalar b42 = -56.0 / 15;
    constexpr Scalar b43 = 32.0 / 9;

    constexpr Scalar b51 = 19372.0 / 6561;
    constexpr Scalar b52 = -25360.0 / 2187;
    constexpr Scalar b53 = 64448.0 / 6561;
    constexpr Scalar b54 = -212.0 / 729;

    constexpr Scalar b61 = 9017.0 / 3168;
    constexpr Scalar b62 = -355.0 / 33;
    constexpr Scalar b63 = 46732.0 / 5247;
    constexpr Scalar b64 = 49.0 / 176;
    constexpr Scalar b65 = -5103.0 / 18656;

    constexpr Scalar c1 = 35.0 / 384;
    // c2 = 0
    constexpr Scalar c3 = 500.0 / 1113;
    constexpr Scalar c4 = 125.0 / 192;
    constexpr Scalar c5 = -2187.0 / 6784;
    constexpr Scalar c6 = 11.0 / 84;

    k1_ = dxdt;  // k1 = 上一次迭代中的 system(x, t)。
    Vector<Scalar, XDim> x = x0 + dt * b21 * k1_;
    k2_ = system.computeFlowMap(t + dt * a2, x);
    x = x0 + dt * b31 * k1_ + dt * b32 * k2_;
    k3_ = system.computeFlowMap(t + dt * a3, x);
    x = x0 + dt * (b41 * k1_ + b42 * k2_ + b43 * k3_);
    k4_ = system.computeFlowMap(t + dt * a4, x);
    x = x0 + dt * (b51 * k1_ + b52 * k2_ + b53 * k3_ + b54 * k4_);
    k5_ = system.computeFlowMap(t + dt * a5, x);
    x = x0 + dt * (b61 * k1_ + b62 * k2_ + b63 * k3_ + b64 * k4_ + b65 * k5_);
    k6_ = system.computeFlowMap(t + dt, x);
    // 更新 x_out 和 dxdt_out（x_out 可以是 x0，dxdt_out 可以是 dxdt）。
    x_out = x0 + dt * (c1 * k1_ + c3 * k3_ + c4 * k4_ + c5 * k5_ + c6 * k6_);
    if (computeDxdtOut) {
      dxdt_out = system.computeFlowMap(t + dt, x_out);
    }
  }

 private:
  /**
   * 估计最大误差值。
   *
   * @param [in] x_old: 上一状态。
   * @param [in] dxdt_old: 上一导数。
   * @param [in] error: 步进误差估计。
   * @param [in] dt: 步长。
   * @param [in] absTol: 绝对误差容忍度。
   * @param [in] relTol: 相对误差容忍度。
   * @return 最大误差值。
   */
  static Scalar maxError(const Vector<Scalar, XDim>& x_old,
                         const Vector<Scalar, XDim>& dxdt_old,
                         const Vector<Scalar, XDim>& x_err, const Scalar dt,
                         const Scalar absTol, const Scalar relTol) {
    const Vector<Scalar, XDim> err =
        x_err.array() /
        (absTol + relTol * (x_old.array().abs() +
                            std::abs(dt) * dxdt_old.array().abs()));
    return err.template lpNorm<Eigen::Infinity>();
  }

  /** Runge-Kutta 步中的中间导数。 */
  Vector<Scalar, XDim> k1_, k2_, k3_, k4_, k5_, k6_;
};

/*
 * 五阶 Runge-Kutta Dormand-Prince（ode45）积分器类。
 *
 * 该实现基于 boost odeint 积分器，并使用
 * controlled boost::numeric::odeint::runge_kutta_dopri5 步进器。
 */
template <typename Scalar, int XDim>
class RungeKuttaDormandPrince5 : public IntegratorBase<Scalar, XDim> {
 public:
  RungeKuttaDormandPrince5() {};

  ~RungeKuttaDormandPrince5() override = default;

  static constexpr size_t maxNumStepsRetries_ = 100;
  /**
   * 基于初始时间、终止时间和步长的等距积分。
   * 长度。
   *
   * @param [in] system: 系统函数。
   * @param [in] observer: 观察器回调。
   * @param [in] initialState: 初始状态。
   * @param [in] startTime: 初始时间。
   * @param [in] finalTime: 终止时间。
   * @param [in] dt: 时间步长。
   */
  void integrateConst(OdeBase<Scalar, XDim>& system,
                      Observer<Scalar, XDim>& observer,
                      const Vector<Scalar, XDim>& initialState,
                      const Scalar startTime, const Scalar finalTime,
                      const Scalar dt) override {
    // 通过加入一个 dt 的小分数确保包含 finalTime，使得：N
    // * dt <= finalTime < (N + 1) * dt。
    const Scalar finalTimeLocal = finalTime + 0.1 * dt;
    const int numSteps = static_cast<int>((finalTimeLocal - startTime) / dt);
    assert(numSteps >= 0);

    Scalar t = startTime;
    Vector<Scalar, XDim> x = initialState;
    Vector<Scalar, XDim> dxdt;
    dxdt = system.computeFlowMap(t, x);
    for (int step = 0; step < numSteps; ++step) {
      observer.observe(x, t);
      const int nextStep = step + 1;
      const Scalar nextTime = startTime + nextStep * dt;
      const bool hasFollowingStep = nextStep < numSteps;
      stepper_.doStep(system, x, dxdt, t, dt, x, dxdt, hasFollowingStep);
      t = nextTime;
    }
    observer.observe(x, t);
  }

 private:
  RungeKuttaDormandPrince5Stepper<Scalar, XDim> stepper_;
};