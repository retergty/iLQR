#include <gtest/gtest.h>

#include <Eigen/LU>
#include <cstdlib>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Tests/Include/MatrixEigenConversion.hpp"
#include "Tests/Include/QpSolver.hpp"
#include "Tests/Include/TestProblemsGeneration.hpp"

using test_tools::matrix_eigen_conversion::fromEigenMatrix;
using test_tools::matrix_eigen_conversion::fromEigenVector;

class QpSolverTest : public testing::Test {
 protected:
  using Scalar = double;
  static constexpr int N_ = 5;
  static constexpr int nx_ = 4;
  static constexpr int nu_ = 3;
  static constexpr int nc_ = 2;
  static constexpr int numDecisionVariables = N_ * (nx_ + nu_) + nx_;
  static constexpr int numConstraints = (N_ + 1) * nx_ + (N_ + 1) * nc_;
  using LinearQuadraticStage_t =
      qp_solver::LinearQuadraticStage<Scalar, nx_, nu_>;
  using CostApproximation_t = qp_solver::QpCostApproximation<Scalar>;
  using ConstraintApproximation_t =
      qp_solver::QpDenseConstraintApproximation<Scalar>;

  QpSolverTest() {
    srand(0);
    lqProblem = generateWellPosedLqProblem();

    cost = qp_solver::getCostMatrices(lqProblem, numDecisionVariables);
    constraints = qp_solver::getConstraintMatrices(
        lqProblem, x0, numConstraints, numDecisionVariables);
    std::tie(primalSolution, dualSolution) =
        qp_solver::solveDenseQp(cost, constraints);
  }

  static LinearQuadraticStage_t generateRandomStage() {
    LinearQuadraticStage_t stage;
    stage.cost.dfdxx = fromEigenMatrix(
        test_tools::getRandomPositiveDefiniteMatrix<Scalar, nx_>());
    stage.cost.dfduu = fromEigenMatrix(
        test_tools::getRandomPositiveDefiniteMatrix<Scalar, nu_>());
    stage.cost.dfdux =
        fromEigenMatrix(Eigen::Matrix<Scalar, nu_, nx_>::Random());
    stage.cost.dfdx = fromEigenVector(Eigen::Vector<Scalar, nx_>::Random());
    stage.cost.dfdu = fromEigenVector(Eigen::Vector<Scalar, nu_>::Random());
    stage.cost.f = std::rand() / static_cast<Scalar>(RAND_MAX);

    stage.dynamics.dfdx =
        fromEigenMatrix(Eigen::Matrix<Scalar, nx_, nx_>::Random());
    stage.dynamics.dfdu =
        fromEigenMatrix(Eigen::Matrix<Scalar, nx_, nu_>::Random());
    stage.dynamics.f = fromEigenVector(Eigen::Vector<Scalar, nx_>::Random());

    stage.constraints.f = Eigen::Vector<Scalar, Eigen::Dynamic>::Random(nc_);
    stage.constraints.dfdx =
        Eigen::Matrix<Scalar, Eigen::Dynamic, nx_>::Random(nc_, nx_);
    stage.constraints.dfdu =
        Eigen::Matrix<Scalar, Eigen::Dynamic, nu_>::Random(nc_, nu_);
    return stage;
  }

  static LinearQuadraticStage_t generateRandomTerminalStage() {
    LinearQuadraticStage_t stage;
    stage.cost.setZero();
    stage.cost.dfdxx = fromEigenMatrix(
        test_tools::getRandomPositiveDefiniteMatrix<Scalar, nx_>());
    stage.cost.dfdx = fromEigenVector(Eigen::Vector<Scalar, nx_>::Random());
    stage.cost.f = std::rand() / static_cast<Scalar>(RAND_MAX);

    stage.dynamics.setZero();
    stage.constraints.f = Eigen::Vector<Scalar, Eigen::Dynamic>::Random(nc_);
    stage.constraints.dfdx =
        Eigen::Matrix<Scalar, Eigen::Dynamic, nx_>::Random(nc_, nx_);
    stage.constraints.dfdu =
        Eigen::Matrix<Scalar, Eigen::Dynamic, nu_>::Random(nc_, nu_);
    return stage;
  }

  static std::vector<LinearQuadraticStage_t> generateRandomLqProblem() {
    std::vector<LinearQuadraticStage_t> lqProblem;
    lqProblem.reserve(N_ + 1);
    for (int k = 0; k < N_; ++k) {
      lqProblem.emplace_back(generateRandomStage());
    }
    lqProblem.emplace_back(generateRandomTerminalStage());
    return lqProblem;
  }

  std::vector<LinearQuadraticStage_t> generateWellPosedLqProblem() {
    for (int retry = 0; retry < 100; ++retry) {
      auto candidate = generateRandomLqProblem();
      const Eigen::Vector<Scalar, nx_> candidateX0 =
          Eigen::Vector<Scalar, nx_>::Random();
      const auto candidateCost =
          qp_solver::getCostMatrices(candidate, numDecisionVariables);
      const auto candidateConstraints = qp_solver::getConstraintMatrices(
          candidate, candidateX0, numConstraints, numDecisionVariables);

      if (candidateCost.dfdxx.fullPivLu().rank() ==
              candidateCost.dfdxx.rows() &&
          candidateConstraints.dfdx.fullPivLu().rank() ==
              candidateConstraints.dfdx.rows()) {
        x0 = candidateX0;
        return candidate;
      }
    }
    throw std::runtime_error(
        "Failed to generate a full-rank random LQ problem.");
  }

  std::vector<LinearQuadraticStage_t> lqProblem;
  Eigen::Vector<Scalar, nx_> x0;
  CostApproximation_t cost;
  ConstraintApproximation_t constraints;
  Eigen::Vector<Scalar, Eigen::Dynamic> primalSolution;
  Eigen::Vector<Scalar, Eigen::Dynamic> dualSolution;
};

TEST_F(QpSolverTest, constraintSatisfaction) {
  ASSERT_TRUE(constraints.f.isApprox(-constraints.dfdx * primalSolution));
}

TEST_F(QpSolverTest, firstOrderOptimality) {
  ASSERT_TRUE(cost.dfdx.isApprox(-cost.dfdxx * primalSolution -
                                 constraints.dfdx.transpose() * dualSolution));
}

TEST_F(QpSolverTest, costMatrixFullRank) {
  ASSERT_TRUE(cost.dfdxx.fullPivLu().rank() == cost.dfdxx.rows());
}

TEST_F(QpSolverTest, constraintMatrixFullRank) {
  ASSERT_TRUE(constraints.dfdx.fullPivLu().rank() == constraints.dfdx.rows());
}
