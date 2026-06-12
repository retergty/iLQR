#include <array>
#include <chrono>
#include <cmath>
#include <iostream>

#include "ExampleModels/ThrustVectorFirstOrderLagConstrain.hpp"
#include "iLQR/iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 30;
constexpr Scalar TimeStep = 0.01;
constexpr Scalar LagAlpha = 0.5;
constexpr size_t LoopCycle = 1000000;

using Descriptor = iLQRDescriptor<
    Scalar, TranscriptionConfig<Dimensions<tvfol_constrain::STATE_DIM,
                                           tvfol_constrain::INPUT_DIM>,
                                Horizon<PredictLength>, DiscreteDynamics>>;
using Solver = iLQR<Descriptor>;
using ThrustVectorILQR =
    tvfol_constrain::ThrustVectorILQR<Scalar, PredictLength>;
using StateVector = typename Solver::StateVector_t;
using InputVector = typename Solver::InputVector_t;

void runOneMpcCycle(ThrustVectorILQR& ilqr,
                    const StateVector& velocityReference, Scalar& currentTime,
                    StateVector& currentState) {
  auto& solver = ilqr.solver();
  auto& problem = ilqr.problem();
  const Vector<Scalar, 3> velocitySetpoint =
      velocityReference.template head<3>();
  const Vector<Scalar, 3> currentVelocity = currentState.template head<3>();
  const Vector<Scalar, 3> commandAccelerationReference =
      currentState.template segment<3>(6);

  ilqr.setDesireTrajectory(currentTime, velocitySetpoint, currentVelocity,
                           commandAccelerationReference);

  solver.run(currentTime, currentState);

  const InputVector firstInput = solver.primalSolution().inputTrajectory_[0];
  currentState = problem.dynamicsPtr->computeMap(currentTime, currentState,
                                                 firstInput, TimeStep);

  currentTime += TimeStep;
}

tvfol_constrain::ThrustVectorILQRSettings<Scalar> makeILQRSettings() {
  tvfol_constrain::ThrustVectorILQRSettings<Scalar> settings;
  settings.ddpSettings.timeStep = TimeStep;
  settings.ddpSettings.maxNumIterations = 20;
  settings.ddpSettings.minRelCost = Scalar(1e-6);

  settings.ocpSettings.Q.setZero();
  settings.ocpSettings.Q.template topLeftCorner<3, 3>().setIdentity();
  settings.ocpSettings.R =
      Scalar(1e-3) * Matrix<Scalar, tvfol_constrain::INPUT_DIM,
                            tvfol_constrain::INPUT_DIM>::Identity();
  settings.ocpSettings.Ra =
      Scalar(1e-4) * Matrix<Scalar, tvfol_constrain::INPUT_DIM,
                            tvfol_constrain::INPUT_DIM>::Identity();
  settings.ocpSettings.Qf.setZero();
  settings.ocpSettings.Qf.template topLeftCorner<3, 3>() =
      Scalar(10.0) * Matrix<Scalar, 3, 3>::Identity();
  settings.ocpSettings.weight = Scalar(1);
  settings.ocpSettings.alpha = LagAlpha;
  return settings;
}
}  // namespace

int main() {
  const auto ilqrSettings = makeILQRSettings();
  auto ilqr = tvfol_constrain::createThrustVectorFirstOrderLagProblem<
      Scalar, PredictLength>(ilqrSettings);

  StateVector currentState;
  currentState.setZero();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(0.0118);
  velocityReference(1) = Scalar(-0.008);
  velocityReference(2) = Scalar(2.45166);

  Scalar currentTime = 0.0;

  std::chrono::time_point now = std::chrono::steady_clock::now();
  for (size_t cycle = 0; cycle < LoopCycle; ++cycle) {
    runOneMpcCycle(ilqr, velocityReference, currentTime, currentState);
  }
  std::cout << "final state: " << currentState.transpose() << std::endl;
  const auto pass_time = std::chrono::steady_clock::now() - now;
  std::cout << "pass time : "
            << std::chrono::duration<double>(pass_time).count() << std::endl;
  std::cout << "average mpc iterations: "
            << ilqr.solver().averageNumIterations() << std::endl;
  return 0;
}
