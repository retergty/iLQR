#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>

#include "LinearSystemDynamics.hpp"
#include "Ocs2QpSolver.hpp"
#include "QuadraticStateCost.hpp"
#include "testProblemsGeneration.hpp"

class Ocs2QpSolverTest : public testing::Test {
protected:
  using Scalar = double;
  static constexpr size_t N = 10;
  static constexpr int STATE_DIM = 2;
  static constexpr int INPUT_DIM = 2;
  static constexpr Scalar precision = 1e-9;
  static constexpr Scalar dt = 1e-2;

  using Problem_t =
      OptimalControlProblem<Scalar, STATE_DIM, INPUT_DIM, N, 0, 0, 0, 0, 0, 0>;
  using Trajectory_t =
      qp_solver::ContinuousTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>;

  void checkDynamics(const Trajectory_t &solution) const {
    Vector<Scalar, STATE_DIM> state = x0;
    for (size_t k = 0; k < N; ++k) {
      const Scalar timeStep =
          solution.timeTrajectory[k + 1] - solution.timeTrajectory[k];
      state +=
          timeStep * system->computeFlowMap(solution.timeTrajectory[k], state,
                                            solution.inputTrajectory[k]);
      ASSERT_TRUE(state.isApprox(solution.stateTrajectory[k + 1], precision));
    }
  }

  Ocs2QpSolverTest() {
    srand(0);

    const auto dynamics =
        test_tools::getRandomDynamics<Scalar, STATE_DIM, INPUT_DIM>();
    A = dynamics.dfdx;
    B = dynamics.dfdu;
    Q = test_tools::getRandomPositiveDefiniteMatrix<Scalar, STATE_DIM>();
    R = test_tools::getRandomPositiveDefiniteMatrix<Scalar, INPUT_DIM>();
    QFinal = test_tools::getRandomPositiveDefiniteMatrix<Scalar, STATE_DIM>();

    system = test_tools::getOcs2Dynamics(dynamics);
    intermediateCost = std::make_unique<
        QuadraticStateInputCost<Scalar, STATE_DIM, INPUT_DIM, N + 1>>(Q, R);
    finalCost =
        std::make_unique<QuadraticStateCost<Scalar, STATE_DIM, N + 1>>(QFinal);

    problem.dynamicsPtr = system.get();
    problem.cost.add(*intermediateCost);
    problem.finalCost.add(*finalCost);

    nominalTrajectory =
        test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
    referenceTrajectory =
        test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
    setReferenceTrajectories(problem, referenceTrajectory);

    x0 = Vector<Scalar, STATE_DIM>::Random();
    solution = qp_solver::solveLinearQuadraticOptimalControlProblem(
        problem, nominalTrajectory, x0);
  }

  void setReferenceTrajectories(Problem_t &targetProblem,
                                const Trajectory_t &trajectory) const {
    targetProblem.timeTrajectory = trajectory.timeTrajectory;
    targetProblem.stateTrajectory = trajectory.stateTrajectory;
    for (size_t k = 0; k < N; ++k) {
      targetProblem.inputTrajectory[k] = trajectory.inputTrajectory[k];
    }
    targetProblem.inputTrajectory[N] = trajectory.inputTrajectory[N - 1];
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
  std::unique_ptr<QuadraticStateCost<Scalar, STATE_DIM, N + 1>> finalCost;
  Problem_t problem;
  Trajectory_t nominalTrajectory;
  Trajectory_t referenceTrajectory;
  Vector<Scalar, STATE_DIM> x0;
  Trajectory_t solution;
};

constexpr size_t Ocs2QpSolverTest::N;
constexpr int Ocs2QpSolverTest::STATE_DIM;
constexpr int Ocs2QpSolverTest::INPUT_DIM;
constexpr Ocs2QpSolverTest::Scalar Ocs2QpSolverTest::precision;
constexpr Ocs2QpSolverTest::Scalar Ocs2QpSolverTest::dt;

TEST_F(Ocs2QpSolverTest, initialCondition) {
  ASSERT_TRUE(x0.isApprox(solution.stateTrajectory.front(), precision));
}

TEST_F(Ocs2QpSolverTest, satisfiesDynamics) { checkDynamics(solution); }

TEST_F(Ocs2QpSolverTest, invariantUnderLinearization) {
  // Different nominalTrajectory, with same time discretization
  auto linearization2 =
      test_tools::getRandomTrajectory<Scalar, STATE_DIM, INPUT_DIM, N>(dt);
  linearization2.timeTrajectory = nominalTrajectory.timeTrajectory;

  // Compare solutions
  auto solution2 = qp_solver::solveLinearQuadraticOptimalControlProblem(
      problem, linearization2, x0);
  for (size_t k = 0; k < N + 1; ++k) {
    ASSERT_TRUE(solution.stateTrajectory[k].isApprox(
        solution2.stateTrajectory[k], precision));
  }
  for (size_t k = 0; k < N; ++k) {
    ASSERT_TRUE(solution.inputTrajectory[k].isApprox(
        solution2.inputTrajectory[k], precision));
  }
}

TEST_F(Ocs2QpSolverTest, knownSolutionAtOrigin) {
  // If the cost's nominal trajectory is set to zero, and the initial state is
  // zero, then the solution has only zeros.
  const auto zeroReference = getZeroTrajectory();
  setReferenceTrajectories(problem, zeroReference);
  const Vector<Scalar, STATE_DIM> zeroX0 = Vector<Scalar, STATE_DIM>::Zero();

  // Obtain solution, with non-zero nominalTrajectory
  auto zeroSolution = qp_solver::solveLinearQuadraticOptimalControlProblem(
      problem, nominalTrajectory, zeroX0);

  for (size_t k = 0; k < N + 1; ++k) {
    ASSERT_TRUE(zeroSolution.stateTrajectory[k].isZero(precision));
  }
  for (size_t k = 0; k < N; ++k) {
    ASSERT_TRUE(zeroSolution.inputTrajectory[k].isZero(precision));
  }
}
