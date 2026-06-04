#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "ExampleModels/CircularKinematics.hpp"
#include "Initialization/DefaultInitializer.hpp"
#include "iLQR/iLQR.hpp"

class CircularKinematicsTest
    : public testing::TestWithParam<circular_model::AugmentedLagrangianType> {
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

  CircularKinematicsTest() {
    for (size_t k = 0; k < PredictLength + 1; ++k) {
      targetTrajectory.timeTrajectory[k] = static_cast<Scalar>(k) * timeStep;
      targetTrajectory.stateTrajectory[k].setZero();
      targetTrajectory.inputTrajectory[k].setZero();
    }
  }

  DDPSettings<Scalar> getSettings() const {
    DDPSettings<Scalar> ddpSettings;
    ddpSettings.timeStep = timeStep;
    ddpSettings.maxNumIterations = 75;
    ddpSettings.minRelCost = 1e-8;
    ddpSettings.lineSearch.minStepLength = 0.01;
    ddpSettings.lineSearch.maxStepLength = 1.0;
    ddpSettings.lineSearch.hessianCorrectionMultiple = 0.01;
    return ddpSettings;
  }

 public:
  static std::string augmentedLagrangianTypeToString(
      circular_model::AugmentedLagrangianType type) {
    switch (type) {
      case circular_model::AugmentedLagrangianType::Quadratic:
        return "QUADRATIC";
      case circular_model::AugmentedLagrangianType::QuadraticStrong:
        return "QUADRATIC_STRONG";
      case circular_model::AugmentedLagrangianType::SmoothAbsolute:
        return "SMOOTH_ABSOLUTE";
      default:
        return "UNKNOWN";
    }
  }

 protected:
  std::string getTestName(
      circular_model::AugmentedLagrangianType augmentedLagrangianType) const {
    return "Circular-Kinematics Test { Algorithm: iLQR, Augmented "
           "Lagrangian: " +
           augmentedLagrangianTypeToString(augmentedLagrangianType) + " }";
  }

  void expectCircularMotion(const Solver_t& solver,
                            const std::string& testName) const {
    const auto performanceIndex = solver.performanceIndex();
    EXPECT_TRUE(std::isfinite(performanceIndex.cost));
    EXPECT_TRUE(std::isfinite(performanceIndex.equalityLagrangian));
    EXPECT_TRUE(std::isfinite(performanceIndex.inequalityLagrangian));
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

    EXPECT_NEAR(solution.stateTrajectory_.back().norm(), Scalar(1.0),
                radiusTolerance)
        << "MESSAGE: " << testName
        << ": failed to keep the trajectory on the unit circle!";
  }

  const Scalar startTime = 0.0;
  const StateVector_t initState{1.0, 0.0};

  TargetTrajectories_t targetTrajectory;
  Initializer_t initializer;
};

TEST_P(CircularKinematicsTest, ILQR) {
  const auto augmentedLagrangianType = GetParam();
  const auto ddpSettings = getSettings();
  auto& problem =
      circular_model::createCircularKinematicsProblem<Scalar, PredictLength>(
          augmentedLagrangianType);

  Solver_t solver(ddpSettings, problem, &initializer);
  solver.setDesireTrajectory(targetTrajectory.timeTrajectory,
                             targetTrajectory.stateTrajectory,
                             targetTrajectory.inputTrajectory);

  ASSERT_NO_THROW(solver.run(startTime, initState));
  expectCircularMotion(solver, getTestName(augmentedLagrangianType));
}

INSTANTIATE_TEST_SUITE_P(
    CircularKinematicsTestCase, CircularKinematicsTest,
    testing::ValuesIn(
        {circular_model::AugmentedLagrangianType::Quadratic,
         circular_model::AugmentedLagrangianType::QuadraticStrong}),
    [](const testing::TestParamInfo<CircularKinematicsTest::ParamType>& info) {
      return CircularKinematicsTest::augmentedLagrangianTypeToString(
          info.param);
    });

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
