#include "Integration.hpp"

/** Runge Kutta Dormand-Prince stepper */
template <typename Scalar, int XDimisions>
class RungeKuttaDormandPrince5Stepper
{
public:
  RungeKuttaDormandPrince5Stepper() = default;
  /**
   * Perform one Dormand-Prince step.
   *
   * @param [in] system: System function.
   * @param [in] x0: current state.
   * @param [in] dxdt: current derivative wrt. time.
   * @param [in] t: current time.
   * @param [in] dt: step size.
   * @param [out] x_out: next state (can be same reference as x0).
   * @param [out] dxdt_out: derivative at next state (can be same reference as dxdt).
   */
  void doStep(const OdeBase<Scalar, XDimisions> &system, const Vector<Scalar, XDimisions> &x0, const Vector<Scalar, XDimisions> &dxdt,
              const Scalar t, const Scalar dt, Vector<Scalar, XDimisions> &x_out, Vector<Scalar, XDimisions> &dxdt_out)
  {
    /* Runge Kutta Dormand-Prince Butcher tableau constants.
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

    k1_ = dxdt; // k1 = system(x, t) from previous iteration
    Vector<Scalar, XDimisions> x = x0 + dt * b21 * k1_;
    k2_ = system.computeFlowMap(t + dt * a2, x);
    x = x0 + dt * b31 * k1_ + dt * b32 * k2_;
    k3_ = system.computeFlowMap(t + dt * a3, x);
    x = x0 + dt * (b41 * k1_ + b42 * k2_ + b43 * k3_);
    k4_ = system.computeFlowMap(t + dt * a4, x);
    x = x0 + dt * (b51 * k1_ + b52 * k2_ + b53 * k3_ + b54 * k4_);
    k5_ = system.computeFlowMap(t + dt * a5, x);
    x = x0 + dt * (b61 * k1_ + b62 * k2_ + b63 * k3_ + b64 * k4_ + b65 * k5_);
    k6_ = system.computeFlowMap(t + dt, x);
    // update x_out and dxdt_out (x_out can be x0 and dxdt_out can be dxdt)
    x_out = x0 + dt * (c1 * k1_ + c3 * k3_ + c4 * k4_ + c5 * k5_ + c6 * k6_);
    dxdt_out = system.computeFlowMap(t + dt, x_out);
  }

private:
  /**
   * Estimate the maximal error value.
   *
   * @param [in] x_old: prevoius state.
   * @param [in] dxdt_old: prevoius derivative.
   * @param [in] error: step error estimate.
   * @param [in] dt: step size.
   * @param [in] absTol: The absolute error tolerance.
   * @param [in] relTol: The relative error tolerance.
   * @return maximal error value.
   */
  static Scalar maxError(const Vector<Scalar, XDimisions> &x_old, const Vector<Scalar, XDimisions> &dxdt_old, const Vector<Scalar, XDimisions> &x_err,
                         const Scalar dt, const Scalar absTol, const Scalar relTol)
  {
    const Vector<Scalar, XDimisions> err = x_err.array() / (absTol + relTol * (x_old.array().abs() + std::abs(dt) * dxdt_old.array().abs()));
    return err.template lpNorm<Eigen::Infinity>();
  }

  /** intermediate derivatives during Runge-Kutta step. */
  Vector<Scalar, XDimisions> k1_, k2_, k3_, k4_, k5_, k6_;
};

/*
 * 5th order Runge Kutta Dormand-Prince (ode45) Integrator class
 *
 * The implementation is based on the boost odeint integrator with the controlled
 * boost::numeric::odeint::runge_kutta_dopri5 stepper.
 */
template <typename Scalar, int XDimisions>
class RungeKuttaDormandPrince5 : public IntegratorBase<Scalar, XDimisions>
{
public:
  RungeKuttaDormandPrince5() {};

  ~RungeKuttaDormandPrince5() override = default;

  static constexpr size_t maxNumStepsRetries_ = 100;
  /**
   * Equidistant integration based on initial and final time as well as step length.
   *
   * @param [in] system: System function
   * @param [in] observer: Observer callback
   * @param [in] initialState: Initial state.
   * @param [in] startTime: Initial time.
   * @param [in] finalTime: Final time.
   * @param [in] dt: Time step.
   */
  void integrateConst(OdeBase<Scalar, XDimisions> &system, Observer<Scalar, XDimisions> &observer,
                      const Vector<Scalar, XDimisions> &initialState, const Scalar startTime, const Scalar finalTime, const Scalar dt) override
  {
    // TODO(mspieler): This does one redundant system() evaluation at the end.

    // Ensure that finalTime is included by adding a fraction of dt such that: N * dt <= finalTime < (N + 1) * dt.
    Scalar finalTimeLocal = finalTime + 0.1 * dt;

    Scalar t = startTime;
    Vector<Scalar, XDimisions> x = initialState;
    Vector<Scalar, XDimisions> dxdt;
    dxdt = system.computeFlowMap(t, x);
    int step = 0;
    while (lessWithSign(t + dt, finalTimeLocal, dt))
    {
      observer.observe(x, t);
      stepper_.doStep(system, x, dxdt, t, dt, x, dxdt);
      step++;
      t = startTime + step * dt;
    }
    observer.observe(x, t);
  }
private:
  /** Helper less comparison for both positive and negative dt case. */
  bool lessWithSign(Scalar t1, Scalar t2, Scalar dt)
  {
    if (dt > 0)
    {
      return t2 - t1 > std::numeric_limits<Scalar>::epsilon();
    }
    else
    {
      return t1 - t2 > std::numeric_limits<Scalar>::epsilon();
    }
  }

private:
  RungeKuttaDormandPrince5Stepper<Scalar, XDimisions> stepper_;
};