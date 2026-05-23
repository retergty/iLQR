#include <gtest/gtest.h>

#include <cmath>

#include "CircularKinematics.hpp"
#include "DefaultInitializer.hpp"
#include "iLQR.hpp"

class CircularKinematicsSmoothAbsolutePenaltyTest : public testing::Test {
 protected:
  using Scalar = double;
  static constexpr int STATE_DIM = circular_model::STATE_DIM;
  static constexpr int INPUT_DIM = circular_model::INPUT_DIM;
  static constexpr Scalar timeStep = 1e-2;
  static constexpr Scalar expectedCost = 0.2;
  static constexpr Scalar radialVelocityTolerance = 1e-1;
  static constexpr Scalar radiusTolerance = 1.2e-1;
  static constexpr Scalar inequalityLagrangianTolerance = 1e-12;
  static constexpr size_t PredictLength = 50;

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

  CircularKinematicsSmoothAbsolutePenaltyTest() {
    for (size_t k = 0; k < PredictLength + 1; ++k) {
      targetTrajectory.timeTrajectory[k] = static_cast<Scalar>(k) * timeStep;
      targetTrajectory.stateTrajectory[k].setZero();
      targetTrajectory.inputTrajectory[k].setZero();
    }
  }

  DDPSettings<Scalar> settings() const {
    DDPSettings<Scalar> ddpSettings;
    ddpSettings.timeStep_ = timeStep;
    ddpSettings.maxNumIterations_ = 40;
    ddpSettings.minRelCost_ = 1e-8;
    ddpSettings.strategy_ = SearchStrategyType::LINE_SEARCH;
    ddpSettings.lineSearch_.minStepLength = 0.01;
    ddpSettings.lineSearch_.maxStepLength = 1.0;
    ddpSettings.lineSearch_.hessianCorrectionMultiple = 1.0;
    return ddpSettings;
  }

  void expectSmoothAbsoluteSolution(const Solver_t& solver) const {
    const auto performanceIndex = solver.performanceIndex();
    EXPECT_TRUE(std::isfinite(performanceIndex.cost));
    EXPECT_TRUE(std::isfinite(performanceIndex.equalityLagrangian));
    EXPECT_TRUE(std::isfinite(performanceIndex.inequalityLagrangian));
    EXPECT_LT(performanceIndex.cost, expectedCost);
    EXPECT_NEAR(performanceIndex.inequalityLagrangian, Scalar(0.0),
                inequalityLagrangianTolerance);

    const auto& solution = solver.primalSolution();
    for (size_t k = 0; k < PredictLength; ++k) {
      const Scalar radialVelocity =
          solution.stateTrajectory_[k].dot(solution.inputTrajectory_[k]);
      EXPECT_NEAR(radialVelocity, Scalar(0.0), radialVelocityTolerance)
          << "k = " << k;
    }

    EXPECT_NEAR(solution.stateTrajectory_.back().norm(), Scalar(1.0),
                radiusTolerance);
  }

  const Scalar startTime = 0.0;
  const StateVector_t initState = (StateVector_t() << 1.0, 0.0).finished();
  TargetTrajectories_t targetTrajectory;
  Initializer_t initializer;
};

TEST_F(CircularKinematicsSmoothAbsolutePenaltyTest,
       ShortHorizonRunsWithSmoothAbsolutePenalty) {
  auto& problem =
      circular_model::createCircularKinematicsProblem<Scalar, PredictLength>(
          circular_model::AugmentedLagrangianType::SmoothAbsolute);
  Solver_t solver(settings(), problem, &initializer);
  solver.setDesireTrajectory(targetTrajectory.timeTrajectory,
                             targetTrajectory.stateTrajectory,
                             targetTrajectory.inputTrajectory);

  ASSERT_NO_THROW(solver.run(startTime, initState));
  expectSmoothAbsoluteSolution(solver);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
