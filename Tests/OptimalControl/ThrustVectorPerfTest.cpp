#include <array>
#include <chrono>
#include <cmath>
#include <iostream>

#include "Initializer.hpp"
#include "ThrustVector.hpp"
#include "iLQR.hpp"

namespace {
using Scalar = double;
constexpr size_t PredictLength = 15;
constexpr Scalar TimeStep = 0.01;
constexpr Scalar Mass = 1.0;
constexpr Scalar Gravity = 9.8;
constexpr size_t LoopCycle = 1000000;

using Descriptor = iLQRDescriptor<
    Scalar, TranscriptionConfig<
                Dimensions<thrust_vector::STATE_DIM, thrust_vector::INPUT_DIM>,
                Horizon<PredictLength>, DiscreteDynamics>>;
using Solver = iLQR<Descriptor>;
using StateVector = typename Solver::StateVector_t;
using InputVector = typename Solver::InputVector_t;
using TimeTrajectory = typename Solver::TimeTrajectory_t;
using StateTrajectory = typename Solver::StateTrajectory_t;
using InputTrajectory = typename Solver::InputTrajectory_t;

class HoverInitializer final
    : public Initializer<Scalar, thrust_vector::STATE_DIM,
                         thrust_vector::INPUT_DIM> {
 public:
  void compute(const Scalar time, const StateVector& state,
               const Scalar nextTime, InputVector& input,
               StateVector& nextState) override {
    (void)time;
    const Scalar dt = nextTime - time;
    input << Scalar(0.0), Scalar(0.0), Mass * Gravity;

    nextState.template head<3>() = state.template head<3>();
    nextState(2) += dt * (input(2) / Mass - Gravity);
    nextState.template tail<3>() = input;
  }
};

void configureVelocityReference(const Scalar initTime,
                                const StateVector& velocityReference,
                                TimeTrajectory& timeTrajectory,
                                StateTrajectory& stateTrajectory,
                                InputTrajectory& inputTrajectory) {
  for (size_t i = 0; i < PredictLength + 1; ++i) {
    timeTrajectory[i] = initTime + static_cast<Scalar>(i) * TimeStep;
    stateTrajectory[i].setZero();
    stateTrajectory[i].template head<3>() =
        velocityReference.template head<3>();
    stateTrajectory[i].template tail<3>() =
        (InputVector() << Scalar(0.0), Scalar(0.0), Mass * Gravity).finished();
    inputTrajectory[i] << Scalar(0.0), Scalar(0.0), Mass * Gravity;
  }
}

void runOneMpcCycle(Solver& solver, Solver::OptimalControlProblem_t& problem,
                    const StateVector& velocityReference, Scalar& currentTime,
                    StateVector& currentState) {
  TimeTrajectory timeTrajectory;
  StateTrajectory stateTrajectory;
  InputTrajectory inputTrajectory;
  configureVelocityReference(currentTime, velocityReference, timeTrajectory,
                             stateTrajectory, inputTrajectory);
  solver.setDesireTrajectory(timeTrajectory, stateTrajectory, inputTrajectory);

  solver.run(currentTime, currentState);
  const InputVector firstInput = solver.primalSolution().inputTrajectory_[0];

  currentState = problem.dynamicsPtr->computeMap(currentTime, currentState,
                                                 firstInput, TimeStep);

  currentTime += TimeStep;
}
}  // namespace

int main() {
  auto& problem =
      thrust_vector::createThrustVectorProblem<Scalar, PredictLength>(Mass);
  HoverInitializer initializer;

  DDPSettings<Scalar> settings;
  settings.timeStep_ = TimeStep;
  settings.maxNumIterations_ = 20;
  settings.minRelCost_ = 1e-6;
  settings.lineSearch_.minStepLength = 0.1;
  Solver solver(settings, problem, &initializer);

  StateVector currentState;
  currentState.setZero();
  currentState.template tail<3>() =
      (InputVector() << Scalar(0.0), Scalar(0.0), Mass * Gravity).finished();

  StateVector velocityReference;
  velocityReference.setZero();
  velocityReference(0) = Scalar(1.0);

  Scalar currentTime = 0.0;

  std::chrono::time_point now = std::chrono::steady_clock::now();
  for (size_t cycle = 0; cycle < LoopCycle; ++cycle) {
    runOneMpcCycle(solver, problem, velocityReference, currentTime,
                   currentState);
  }
  std::cout << "final state: " << currentState.transpose() << std::endl;
  const auto pass_time = std::chrono::steady_clock::now() - now;
  std::cout << "pass time : "
            << std::chrono::duration<double>(pass_time).count() << std::endl;
  return 0;
}
