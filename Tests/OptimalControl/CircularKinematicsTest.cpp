#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "CircularKinematics.hpp"
#include "DefaultInitializer.hpp"
#include "iLQR.hpp"

class CircularKinematicsTest
    : public testing::TestWithParam<SearchStrategyType> {
 protected:
  using Scalar = double;
  static constexpr int STATE_DIM = circular_model::STATE_DIM;
  static constexpr int INPUT_DIM = circular_model::INPUT_DIM;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar expectedCost = 0.1;
  static constexpr Scalar radialVelocityTolerance = 2e-2;
  static constexpr Scalar radiusTolerance = 5e-2;
  static constexpr Scalar equalityLagrangianTolerance = 3e-2;
  static constexpr Scalar inequalityLagrangianTolerance = 1e-12;
  static constexpr size_t PredictLength = 1000;

  using ConstraintConfig_t = ConstraintConfig<
      StateConstraintConfig<ConstraintLayout<>>,
      StateInputConstraintConfig<ConstraintLayout<
          ConstraintGroupLayout<ConstraintTerm<1>>, ConstraintGroupLayout<>>>,
      FinalStateConstraintConfig<ConstraintLayout<>>>;
  using Descriptor_t =
      iLQRDescriptor<Scalar,
                     TranscriptionConfig<Dimensions<STATE_DIM, INPUT_DIM>,
                                         Horizon<PredictLength>>,
                     ConstraintConfig_t>;
  using Solver_t = iLQR<Descriptor_t>;
  using Initializer_t = DefaultInitializer<Scalar, STATE_DIM, INPUT_DIM>;
  using StateVector_t = typename Solver_t::StateVector_t;
  using TargetTrajectories_t = typename Solver_t::TargetTrajectories_t;

  CircularKinematicsTest()
      : problem(
            circular_model::createCircularKinematicsProblem<Scalar,
                                                            PredictLength>()) {
    setZeroReferenceTrajectory();
  }

  DDPSettings<Scalar> getSettings(SearchStrategyType strategy) const {
    DDPSettings<Scalar> ddpSettings;
    ddpSettings.timeStep_ = timeStep;
    ddpSettings.maxNumIterations_ = 75;
    ddpSettings.minRelCost_ = 1e-8;
    ddpSettings.strategy_ = strategy;
    ddpSettings.lineSearch_.minStepLength = 0.01;
    ddpSettings.lineSearch_.maxStepLength = 1.0;
    ddpSettings.lineSearch_.hessianCorrectionMultiple = 0.01;
    return ddpSettings;
  }

  std::string getTestName(SearchStrategyType strategy) const {
    switch (strategy) {
      case SearchStrategyType::LINE_SEARCH:
        return "Circular-Kinematics Test { Algorithm: iLQR, Strategy: "
               "LINE_SEARCH }";
      case SearchStrategyType::LEVENBERG_MARQUARDT:
        return "Circular-Kinematics Test { Algorithm: iLQR, Strategy: "
               "LEVENBERG_MARQUARDT }";
      default:
        return "Circular-Kinematics Test { Algorithm: iLQR, Strategy: "
               "UNKNOWN }";
    }
  }

  void setZeroReferenceTrajectory() {
    for (size_t k = 0; k < PredictLength + 1; ++k) {
      targetTrajectory.timeTrajectory[k] = static_cast<Scalar>(k) * timeStep;
      targetTrajectory.stateTrajectory[k].setZero();
      targetTrajectory.inputTrajectory[k].setZero();
    }
  }

  void expectCircularMotion(const Solver_t& solver,
                            const std::string& testName) const {
    const auto performanceIndex = solver.performanceIndex();
    EXPECT_LT(performanceIndex.cost, expectedCost)
        << "MESSAGE: " << testName << ": failed in the cost test!";
    EXPECT_NEAR(performanceIndex.equalityLagrangian, Scalar(0.0),
                equalityLagrangianTolerance)
        << "MESSAGE: " << testName
        << ": failed in the equality Lagrangian test!";
    EXPECT_NEAR(performanceIndex.inequalityLagrangian, Scalar(0.0),
                inequalityLagrangianTolerance)
        << "MESSAGE: " << testName
        << ": failed in the inequality Lagrangian test!";

    const auto& solution = solver.primalSolution();
    for (size_t k = 0; k < PredictLength; ++k) {
      const Scalar radialVelocity =
          solution.stateTrajectory_[k].dot(solution.inputTrajectory_[k]);
      EXPECT_NEAR(radialVelocity, Scalar(0.0), radialVelocityTolerance)
          << "MESSAGE: " << testName
          << ": failed in state-input equality constraint at k = " << k;
    }

    const Scalar finalRadius = solution.stateTrajectory_.back().norm();
    EXPECT_NEAR(finalRadius, Scalar(1.0), radiusTolerance)
        << "MESSAGE: " << testName
        << ": failed to keep the trajectory on the unit circle!";
  }

  const Scalar startTime = 0.0;
  const StateVector_t initState = (StateVector_t() << 1.0, 0.0).finished();

  typename Solver_t::OptimalControlProblem_t& problem;
  TargetTrajectories_t targetTrajectory;
  Initializer_t initializer;
};

constexpr int CircularKinematicsTest::STATE_DIM;
constexpr int CircularKinematicsTest::INPUT_DIM;
constexpr CircularKinematicsTest::Scalar CircularKinematicsTest::timeStep;
constexpr CircularKinematicsTest::Scalar CircularKinematicsTest::expectedCost;
constexpr CircularKinematicsTest::Scalar
    CircularKinematicsTest::radialVelocityTolerance;
constexpr CircularKinematicsTest::Scalar
    CircularKinematicsTest::radiusTolerance;
constexpr CircularKinematicsTest::Scalar
    CircularKinematicsTest::equalityLagrangianTolerance;
constexpr CircularKinematicsTest::Scalar
    CircularKinematicsTest::inequalityLagrangianTolerance;
constexpr size_t CircularKinematicsTest::PredictLength;

TEST_P(CircularKinematicsTest, ILQR) {
  const auto strategy = GetParam();
  const auto ddpSettings = getSettings(strategy);

  Solver_t solver(ddpSettings, problem, &initializer);
  solver.setDesireTrajectory(targetTrajectory.timeTrajectory,
                             targetTrajectory.stateTrajectory,
                             targetTrajectory.inputTrajectory);

  ASSERT_NO_THROW(solver.run(startTime, initState));
  expectCircularMotion(solver, getTestName(strategy));
}

INSTANTIATE_TEST_SUITE_P(
    CircularKinematicsTestCase, CircularKinematicsTest,
    testing::ValuesIn({SearchStrategyType::LINE_SEARCH}),
    [](const testing::TestParamInfo<CircularKinematicsTest::ParamType>& info) {
      switch (info.param) {
        case SearchStrategyType::LINE_SEARCH:
          return std::string("LINE_SEARCH");
        case SearchStrategyType::LEVENBERG_MARQUARDT:
          return std::string("LEVENBERG_MARQUARDT");
        default:
          return std::string("UNKNOWN");
      }
    });

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
