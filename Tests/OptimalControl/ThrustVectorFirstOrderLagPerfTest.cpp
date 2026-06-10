#include <array>
#include <chrono>
#include <cmath>
#include <iostream>

#include "ExampleModels/ThrustVectorFirstOrderLag.hpp"
#include "iLQR/iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 30;
constexpr Scalar TimeStep = 0.01;
constexpr Scalar LagAlpha = 0.5;
constexpr size_t LoopCycle = 1000000;

using Descriptor = iLQRDescriptor<
    Scalar,
    TranscriptionConfig<Dimensions<thrust_vector_first_order_lag::STATE_DIM,
                                   thrust_vector_first_order_lag::INPUT_DIM>,
                        Horizon<PredictLength>, DiscreteDynamics>>;
using Solver = iLQR<Descriptor>;
using ThrustVectorILQR =
    thrust_vector_first_order_lag::ThrustVectorILQR<Scalar, PredictLength>;
using StateVector = typename Solver::StateVector_t;
using InputVector = typename Solver::InputVector_t;

size_t runOneMpcCycle(ThrustVectorILQR& ilqr,
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

  const size_t iterationsBefore = solver.totalNumIterations();
  solver.run(currentTime, currentState);
  const size_t mpcIterations = solver.totalNumIterations() - iterationsBefore;

  const InputVector firstInput = solver.primalSolution().inputTrajectory_[0];
  currentState = problem.dynamicsPtr->computeMap(currentTime, currentState,
                                                 firstInput, TimeStep);

  currentTime += TimeStep;
  return mpcIterations;
}

thrust_vector_first_order_lag::ThrustVectorILQRSettings<Scalar>
makeILQRSettings() {
  thrust_vector_first_order_lag::ThrustVectorILQRSettings<Scalar> settings;
  settings.ddpSettings.timeStep = TimeStep;
  settings.ddpSettings.maxNumIterations = 20;
  settings.ddpSettings.minRelCost = Scalar(1e-6);

  settings.ocpSettings.Q.setZero();
  settings.ocpSettings.Q.template topLeftCorner<3, 3>().setIdentity();
  settings.ocpSettings.R =
      Scalar(1e-3) *
      Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
             thrust_vector_first_order_lag::INPUT_DIM>::Identity();
  settings.ocpSettings.Ra =
      Scalar(1e-4) *
      Matrix<Scalar, thrust_vector_first_order_lag::INPUT_DIM,
             thrust_vector_first_order_lag::INPUT_DIM>::Identity();
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
  auto ilqr =
      thrust_vector_first_order_lag::createThrustVectorFirstOrderLagProblem<
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
  size_t totalMpcIterations = 0;
  for (size_t cycle = 0; cycle < LoopCycle; ++cycle) {
    totalMpcIterations +=
        runOneMpcCycle(ilqr, velocityReference, currentTime, currentState);
  }
  std::cout << "final state: " << currentState.transpose() << std::endl;
  const auto pass_time = std::chrono::steady_clock::now() - now;
  std::cout << "pass time : "
            << std::chrono::duration<double>(pass_time).count() << std::endl;
  std::cout << "average mpc iterations: "
            << static_cast<double>(totalMpcIterations) /
                   static_cast<double>(LoopCycle)
            << std::endl;
  return 0;
}
