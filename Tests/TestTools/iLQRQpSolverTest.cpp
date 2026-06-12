#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>

#include "Cost/QuadraticStateCost.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "Integration/SensitivityIntegrator.hpp"
#include "Tests/Include/MatrixEigenConversion.hpp"
#include "Tests/Include/TestProblemsGeneration.hpp"
#include "Tests/Include/iLQRQpSolver.hpp"

using test_tools::matrix_eigen_conversion::fromEigenMatrix;
using test_tools::matrix_eigen_conversion::fromEigenVector;
using test_tools::matrix_eigen_conversion::toEigenMatrix;
using test_tools::matrix_eigen_conversion::toEigenVector;

class iLQRQpSolverTest : public testing::Test {
 protected:
  using Scalar = double;
  static constexpr size_t N = 10;
  static constexpr int STATE_DIM = 2;
  static constexpr int INPUT_DIM = 2;
  static constexpr Scalar precision = 1e-9;
  static constexpr Scalar dt = 1e-2;

  using Problem_t =
      qp_solver::QpContinuousDynamicsOptimalControlProblem_t<Scalar, STATE_DIM,
                                                             INPUT_DIM, N>;
  using Trajectory_t =
      qp_solver::ContinuousTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>;
  using Transcription_t =
      qp_solver::QpContinuousDynamicsTranscriptionConfig_t<Scalar, STATE_DIM,
                                                           INPUT_DIM, N>;
  using TargetTrajectories_t = TargetTrajectories<Scalar, Transcription_t>;

  void checkDeviationDynamics(const Trajectory_t& solution,
                              const Trajectory_t& nominal) {
    const auto lqp = qp_solver::getLinearQuadraticApproximation(
        problem, targetTrajectory, nominal);

    ASSERT_TRUE(
        (solution.stateTrajectory.front() - nominal.stateTrajectory.front())
            .isApprox(x0 - nominal.stateTrajectory.front(), precision));

    for (size_t k = 0; k < N; ++k) {
      const auto dx = solution.stateTrajectory[k] - nominal.stateTrajectory[k];
      const auto du = solution.inputTrajectory[k] - nominal.inputTrajectory[k];
      const Eigen::Vector<Scalar, STATE_DIM> expectedNextDx =
          toEigenMatrix(lqp[k].dynamics.dfdx) * dx +
          toEigenMatrix(lqp[k].dynamics.dfdu) * du +
          toEigenVector(lqp[k].dynamics.f);
      const auto nextDx =
          solution.stateTrajectory[k + 1] - nominal.stateTrajectory[k + 1];
      ASSERT_TRUE(expectedNextDx.isApprox(nextDx, precision));
    }
  }

  Trajectory_t getDynamicallyConsistentTrajectory() const {
    EK2DynamicsDiscretizer<Scalar, STATE_DIM, INPUT_DIM> discretizer;
    Trajectory_t trajectory;

    trajectory.stateTrajectory[0] = Eigen::Vector<Scalar, STATE_DIM>::Random();
    for (size_t k = 0; k < N; ++k) {
      trajectory.timeTrajectory[k] = static_cast<Scalar>(k) * dt;
      trajectory.inputTrajectory[k] =
          Eigen::Vector<Scalar, INPUT_DIM>::Random();
      const auto nextState = discretizer.discretize(
          *system, trajectory.timeTrajectory[k],
          fromEigenVector(trajectory.stateTrajectory[k]),
          fromEigenVector(trajectory.inputTrajectory[k]), dt);
      trajectory.stateTrajectory[k + 1] = toEigenVector(nextState);
    }
    trajectory.timeTrajectory[N] = static_cast<Scalar>(N) * dt;

    return trajectory;
  }

  iLQRQpSolverTest() {
    srand(0);

    const auto dynamics =
        test_tools::getRandomDynamics<Scalar, STATE_DIM, INPUT_DIM>();
    A = dynamics.dfdx;
    B = dynamics.dfdu;
    Q = fromEigenMatrix(
        test_tools::getRandomPositiveDefiniteMatrix<Scalar, STATE_DIM>());
    R = fromEigenMatrix(
        test_tools::getRandomPositiveDefiniteMatrix<Scalar, INPUT_DIM>());
    QFinal = fromEigenMatrix(
        test_tools::getRandomPositiveDefiniteMatrix<Scalar, STATE_DIM>());

    system = test_tools::getiLQRDynamics(dynamics);
    intermediateCost = std::make_unique<
        QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>(Q, R, 0);
    finalCost = std::make_unique<
        QuadraticStateCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>(QFinal, 0);

    problem.dynamicsPtr = system.get();
    problem.cost.add(*intermediateCost);
    problem.finalCost.add(*finalCost);

    nominalTrajectory =
        test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
    referenceTrajectory =
        test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
    setReferenceTrajectories(referenceTrajectory);

    x0 = Eigen::Vector<Scalar, STATE_DIM>::Random();
    solution = qp_solver::solveLinearQuadraticOptimalControlProblem(
        problem, targetTrajectory, nominalTrajectory, x0);
  }

  void setReferenceTrajectories(const Trajectory_t& trajectory) {
    targetTrajectory.timeTrajectory = trajectory.timeTrajectory;
    for (size_t k = 0; k < N + 1; ++k) {
      targetTrajectory.stateTrajectory[k] =
          fromEigenVector(trajectory.stateTrajectory[k]);
    }
    for (size_t k = 0; k < N; ++k) {
      targetTrajectory.inputTrajectory[k] =
          fromEigenVector(trajectory.inputTrajectory[k]);
    }
    targetTrajectory.inputTrajectory[N] =
        fromEigenVector(trajectory.inputTrajectory[N - 1]);
  }

  static Trajectory_t getZeroTrajectory() {
    Trajectory_t trajectory;
    for (size_t k = 0; k < N + 1; ++k) {
      trajectory.timeTrajectory[k] = static_cast<Scalar>(k) * dt;
      trajectory.stateTrajectory[k].setZero();
    }
    for (size_t k = 0; k < N; ++k) {
      trajectory.inputTrajectory[k].setZero();
    }
    return trajectory;
  }

  Matrix<Scalar, STATE_DIM, STATE_DIM> A;
  Matrix<Scalar, STATE_DIM, INPUT_DIM> B;
  Matrix<Scalar, STATE_DIM, STATE_DIM> Q;
  Matrix<Scalar, INPUT_DIM, INPUT_DIM> R;
  Matrix<Scalar, STATE_DIM, STATE_DIM> QFinal;

  std::unique_ptr<LinearSystemDynamics<Scalar, STATE_DIM, INPUT_DIM>> system;
  std::unique_ptr<QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>
      intermediateCost;
  std::unique_ptr<QuadraticStateCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>
      finalCost;
  Problem_t problem;
  TargetTrajectories_t targetTrajectory;
  Trajectory_t nominalTrajectory;
  Trajectory_t referenceTrajectory;
  Eigen::Vector<Scalar, STATE_DIM> x0;
  Trajectory_t solution;
};

TEST_F(iLQRQpSolverTest, initialCondition) {
  ASSERT_TRUE(x0.isApprox(solution.stateTrajectory.front(), precision));
}

TEST_F(iLQRQpSolverTest, satisfiesDeviationDynamics) {
  checkDeviationDynamics(solution, nominalTrajectory);
}

TEST_F(iLQRQpSolverTest, invariantUnderLinearization) {
  // 对于具有二次代价的线性系统，绝对解对动态一致的
  // 名义轨迹保持不变。
  const auto linearization1 = getDynamicallyConsistentTrajectory();
  const auto linearization2 = getDynamicallyConsistentTrajectory();

  // 比较解。
  const auto solution1 = qp_solver::solveLinearQuadraticOptimalControlProblem(
      problem, targetTrajectory, linearization1, x0);
  const auto solution2 = qp_solver::solveLinearQuadraticOptimalControlProblem(
      problem, targetTrajectory, linearization2, x0);
  for (size_t k = 0; k < N + 1; ++k) {
    ASSERT_TRUE(solution1.stateTrajectory[k].isApprox(
        solution2.stateTrajectory[k], precision));
  }
  for (size_t k = 0; k < N; ++k) {
    ASSERT_TRUE(solution1.inputTrajectory[k].isApprox(
        solution2.inputTrajectory[k], precision));
  }
}

TEST_F(iLQRQpSolverTest, knownSolutionAtOrigin) {
  // 如果代价的名义轨迹设为零且初始状态为
  // 零，则解也全为零。
  const auto zeroReference = getZeroTrajectory();
  setReferenceTrajectories(zeroReference);
  const Eigen::Vector<Scalar, STATE_DIM> zeroX0 =
      Eigen::Vector<Scalar, STATE_DIM>::Zero();

  // 获取零名义轨迹下的解。
  auto zeroSolution = qp_solver::solveLinearQuadraticOptimalControlProblem(
      problem, targetTrajectory, zeroReference, zeroX0);

  for (size_t k = 0; k < N + 1; ++k) {
    ASSERT_TRUE(zeroSolution.stateTrajectory[k].isZero(precision));
  }
  for (size_t k = 0; k < N; ++k) {
    ASSERT_TRUE(zeroSolution.inputTrajectory[k].isZero(precision));
  }
}
