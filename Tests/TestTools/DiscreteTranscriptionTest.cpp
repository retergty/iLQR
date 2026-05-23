
#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>
#include <vector>

#include "LinearSystemDynamics.hpp"
#include "QpDiscreteTranscription.hpp"
#include "QuadraticStateCost.hpp"
#include "testProblemsGeneration.hpp"

class DiscreteTranscriptionTest : public testing::Test {
 protected:
  using Scalar = double;
  static constexpr size_t N = 10;
  static constexpr int STATE_DIM = 3;
  static constexpr int INPUT_DIM = 2;
  static constexpr Scalar dt = 1e-2;

  using Problem_t =
      qp_solver::QpOptimalControlProblem_t<Scalar, STATE_DIM, INPUT_DIM, N, 0,
                                           0, 0, 0, 0, 0>;
  using Trajectory_t =
      qp_solver::ContinuousTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>;
  using Transcription_t =
      qp_solver::QpTranscriptionConfig_t<Scalar, STATE_DIM, INPUT_DIM, N, 0, 0,
                                         0, 0, 0, 0>;
  using TargetTrajectories_t = TargetTrajectories<Scalar, Transcription_t>;
  using LinearQuadraticStage_t =
      qp_solver::LinearQuadraticStage<Scalar, STATE_DIM, INPUT_DIM>;

  DiscreteTranscriptionTest() {
    srand(0);

    A << 1.0, 0.1, -0.2, 0.0, -0.5, 0.3, 0.4, 0.0, 0.2;
    B << 1.0, 0.0, 0.0, 1.0, 0.5, -0.5;
    Q << 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0;
    QState << 0.5, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.75;
    R << 5.0, 0.0, 0.0, 6.0;
    QFinal << 7.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 9.0;

    system =
        std::make_unique<LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>>(A,
                                                                             B);
    intermediateCost = std::make_unique<
        QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>(Q, R);
    stateCost =
        std::make_unique<QuadraticStateCost<Scalar, STATE_DIM, N + 1>>(QState);
    finalCost =
        std::make_unique<QuadraticStateCost<Scalar, STATE_DIM, N + 1>>(QFinal);

    problem.dynamicsPtr = system.get();
    problem.cost.add(*intermediateCost);
    problem.stateCost.add(*stateCost);
    problem.finalCost.add(*finalCost);

    linearization =
        test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
    setReferenceTrajectories(linearization);

    unconstrainedLqr = qp_solver::getLinearQuadraticApproximation(
        problem, targetTrajectory, linearization);
  }

  void setReferenceTrajectories(const Trajectory_t& trajectory) {
    targetTrajectory.timeTrajectory = trajectory.timeTrajectory;
    targetTrajectory.stateTrajectory = trajectory.stateTrajectory;
    for (size_t k = 0; k < N; ++k) {
      targetTrajectory.inputTrajectory[k] = trajectory.inputTrajectory[k];
    }
    targetTrajectory.inputTrajectory[N] = trajectory.inputTrajectory[N - 1];
  }

  void checkSizes(const std::vector<LinearQuadraticStage_t>& lqr) const {
    ASSERT_EQ(lqr.size(), N + 1);
    for (size_t k = 0; k < N; ++k) {
      // Cost sizes
      ASSERT_EQ(lqr[k].cost.dfdxx.rows(), STATE_DIM);
      ASSERT_EQ(lqr[k].cost.dfdxx.cols(), STATE_DIM);
      ASSERT_EQ(lqr[k].cost.dfdux.rows(), INPUT_DIM);
      ASSERT_EQ(lqr[k].cost.dfdux.cols(), STATE_DIM);
      ASSERT_EQ(lqr[k].cost.dfduu.rows(), INPUT_DIM);
      ASSERT_EQ(lqr[k].cost.dfduu.cols(), INPUT_DIM);

      // Dynamics sizes
      ASSERT_EQ(lqr[k].dynamics.dfdx.rows(), STATE_DIM);
      ASSERT_EQ(lqr[k].dynamics.dfdx.cols(), STATE_DIM);
      ASSERT_EQ(lqr[k].dynamics.dfdu.rows(), STATE_DIM);
      ASSERT_EQ(lqr[k].dynamics.dfdu.cols(), INPUT_DIM);

      // Constraint sizes
      ASSERT_EQ(lqr[k].constraints.f.rows(), 0);
      ASSERT_EQ(lqr[k].constraints.dfdx.rows(), 0);
      ASSERT_EQ(lqr[k].constraints.dfdx.cols(), STATE_DIM);
      ASSERT_EQ(lqr[k].constraints.dfdu.rows(), 0);
      ASSERT_EQ(lqr[k].constraints.dfdu.cols(), INPUT_DIM);
    }

    // Terminal Cost size
    ASSERT_EQ(lqr[N].cost.dfdxx.rows(), STATE_DIM);
    ASSERT_EQ(lqr[N].cost.dfdxx.cols(), STATE_DIM);

    // Terminal Constraint size
    ASSERT_EQ(lqr[N].constraints.f.rows(), 0);
    ASSERT_EQ(lqr[N].constraints.dfdx.rows(), 0);
    ASSERT_EQ(lqr[N].constraints.dfdx.cols(), STATE_DIM);
  }

  Matrix<Scalar, STATE_DIM, STATE_DIM> A;
  Matrix<Scalar, STATE_DIM, INPUT_DIM> B;
  Matrix<Scalar, STATE_DIM, STATE_DIM> Q;
  Matrix<Scalar, STATE_DIM, STATE_DIM> QState;
  Matrix<Scalar, INPUT_DIM, INPUT_DIM> R;
  Matrix<Scalar, STATE_DIM, STATE_DIM> QFinal;

  std::unique_ptr<LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>> system;
  std::unique_ptr<QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>
      intermediateCost;
  std::unique_ptr<QuadraticStateCost<Scalar, STATE_DIM, N + 1>> stateCost;
  std::unique_ptr<QuadraticStateCost<Scalar, STATE_DIM, N + 1>> finalCost;
  Problem_t problem;
  TargetTrajectories_t targetTrajectory;
  Trajectory_t linearization;
  std::vector<LinearQuadraticStage_t> unconstrainedLqr;
};

TEST_F(DiscreteTranscriptionTest, unconstrainedLqrHasCorrectSizes) {
  checkSizes(unconstrainedLqr);
}

TEST_F(DiscreteTranscriptionTest, ek2DeviationDynamicsMatchesLinearSystem) {
  const auto expectedA = Matrix<Scalar, STATE_DIM, STATE_DIM>::Identity() +
                         A * dt + Scalar(0.5) * A * A * dt * dt;
  const auto expectedB = B * dt + Scalar(0.5) * A * B * dt * dt;

  for (size_t k = 0; k < N; ++k) {
    ASSERT_TRUE(unconstrainedLqr[k].dynamics.dfdx.isApprox(expectedA));
    ASSERT_TRUE(unconstrainedLqr[k].dynamics.dfdu.isApprox(expectedB));
    ASSERT_TRUE(unconstrainedLqr[k].dynamics.f.isZero());
  }
}

TEST_F(DiscreteTranscriptionTest, costHessiansHaveExpectedScaling) {
  const auto expectedQ = (Q + QState) * dt;
  const auto expectedR = R * dt;

  for (size_t k = 0; k < N; ++k) {
    ASSERT_TRUE(unconstrainedLqr[k].cost.dfdxx.isApprox(expectedQ));
    ASSERT_TRUE(unconstrainedLqr[k].cost.dfduu.isApprox(expectedR));
    ASSERT_TRUE(unconstrainedLqr[k].cost.dfdux.isZero());
  }

  ASSERT_TRUE(unconstrainedLqr[N].cost.dfdxx.isApprox(QFinal));
}

TEST_F(DiscreteTranscriptionTest, linearizationInvariance) {
  auto linearization2 =
      test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
  linearization2.timeTrajectory = linearization.timeTrajectory;

  setReferenceTrajectories(linearization2);
  const auto lqp2 = qp_solver::getLinearQuadraticApproximation(
      problem, targetTrajectory, linearization2);

  // All Hessians and Jacobians stay the same. The linear and constant parts
  // change with the nominal trajectory.
  for (size_t k = 0; k < N; ++k) {
    // Cost
    ASSERT_TRUE(unconstrainedLqr[k].cost.dfdxx.isApprox(lqp2[k].cost.dfdxx));
    ASSERT_TRUE(unconstrainedLqr[k].cost.dfdux.isApprox(lqp2[k].cost.dfdux));
    ASSERT_TRUE(unconstrainedLqr[k].cost.dfduu.isApprox(lqp2[k].cost.dfduu));

    // Dynamics
    ASSERT_TRUE(
        unconstrainedLqr[k].dynamics.dfdx.isApprox(lqp2[k].dynamics.dfdx));
    ASSERT_TRUE(
        unconstrainedLqr[k].dynamics.dfdu.isApprox(lqp2[k].dynamics.dfdu));

    // Constraints
    ASSERT_TRUE(unconstrainedLqr[k].constraints.dfdx.isApprox(
        lqp2[k].constraints.dfdx));
    ASSERT_TRUE(unconstrainedLqr[k].constraints.dfdu.isApprox(
        lqp2[k].constraints.dfdu));
  }

  // Terminal Cost
  ASSERT_TRUE(unconstrainedLqr[N].cost.dfdxx.isApprox(lqp2[N].cost.dfdxx));

  // Terminal Constraints
  ASSERT_TRUE(
      unconstrainedLqr[N].constraints.dfdx.isApprox(lqp2[N].constraints.dfdx));
}
