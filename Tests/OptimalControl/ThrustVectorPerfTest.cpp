#include <array>
#include <chrono>
#include <cmath>
#include <iostream>

#include "ExampleModels/ThrustVector.hpp"
#include "Initialization/Initializer.hpp"
#include "iLQR/iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 30;
constexpr Scalar TimeStep = 0.01;
constexpr size_t LoopCycle = 1000000;

using Descriptor = iLQRDescriptor<
    Scalar, TranscriptionConfig<
                Dimensions<thrust_vector::STATE_DIM, thrust_vector::INPUT_DIM>,
                Horizon<PredictLength>, DiscreteDynamics>>;
using Solver = iLQR<Descriptor>;
using StateVector = typename Solver::StateVector_t;
using InputVector = typename Solver::InputVector_t;
using TargetTrajectories = typename Solver::TargetTrajectories_t;

InputVector hoverInput() {
  return InputVector{Scalar(0.0), Scalar(0.0), Scalar(0.0)};
}

class HoverInitializer final
    : public Initializer<Scalar, thrust_vector::STATE_DIM,
                         thrust_vector::INPUT_DIM> {
 public:
  void compute(const Scalar time, const StateVector& state,
               const Scalar nextTime, InputVector& input,
               StateVector& nextState) override {
    (void)time;
    const Scalar dt = nextTime - time;
    input = hoverInput();

    const InputVector currentAcceleration =
        state.template tail<thrust_vector::INPUT_DIM>() + input;
    nextState.template head<3>() =
        state.template head<3>() + dt * currentAcceleration;
    nextState.template tail<3>() = currentAcceleration;
  }
};

void configureVelocityReference(const Scalar initTime,
                                const StateVector& velocityReference,
                                TargetTrajectories& targetTrajectory) {
  for (size_t i = 0; i < PredictLength + 1; ++i) {
    targetTrajectory.timeTrajectory[i] =
        initTime + static_cast<Scalar>(i) * TimeStep;
    targetTrajectory.stateTrajectory[i].setZero();
    targetTrajectory.stateTrajectory[i].template head<3>() =
        velocityReference.template head<3>();
    targetTrajectory.stateTrajectory[i].template tail<3>() = hoverInput();
    targetTrajectory.inputTrajectory[i] = hoverInput();
  }
}

size_t runOneMpcCycle(Solver& solver, Solver::OptimalControlProblem_t& problem,
                      const StateVector& velocityReference, Scalar& currentTime,
                      StateVector& currentState) {
  configureVelocityReference(currentTime, velocityReference,
                             solver.targetTrajectory());

  const size_t iterationsBefore = solver.totalNumIterations();
  solver.run(currentTime, currentState);
  const size_t mpcIterations = solver.totalNumIterations() - iterationsBefore;

  const InputVector firstInput = solver.primalSolution().inputTrajectory_[0];

  currentState = problem.dynamicsPtr->computeMap(currentTime, currentState,
                                                 firstInput, TimeStep);

  currentTime += TimeStep;
  return mpcIterations;
}
}  // namespace

int main() {
  auto& problem =
      thrust_vector::createThrustVectorProblem<Scalar, PredictLength>();
  HoverInitializer initializer;

  DDPSettings<Scalar> settings;
  settings.timeStep = TimeStep;
  settings.maxNumIterations = 20;
  settings.minRelCost = 1e-6;
  Solver solver(settings, problem, &initializer);

  StateVector currentState;
  currentState.setZero();
  currentState.template tail<3>() = hoverInput();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(0.0118);
  velocityReference(1) = Scalar(-0.008);
  velocityReference(2) = Scalar(2.45166);

  Scalar currentTime = 0.0;

  std::chrono::time_point now = std::chrono::steady_clock::now();
  size_t totalMpcIterations = 0;
  for (size_t cycle = 0; cycle < LoopCycle; ++cycle) {
    totalMpcIterations += runOneMpcCycle(solver, problem, velocityReference,
                                         currentTime, currentState);
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
