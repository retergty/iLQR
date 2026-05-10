/******************************************************************************
Copyright (c) 2017, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

#pragma once

#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"
#include <Eigen/LU>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace qp_solver {

  template <typename Scalar, int XDim, int UDim>
  VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0> getConstraintMatrices(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
    const Vector<Scalar, XDim>& dx0, int numConstraints, int numDecisionVariables);

  template <typename Scalar, int XDim, int UDim>
  ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0> getCostMatrices(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp, int numDecisionVariables);

  template <typename Scalar>
  std::pair<Vector<Scalar, Eigen::Dynamic>, Vector<Scalar, Eigen::Dynamic>> solveDenseQp(
    const ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>& cost,
    const VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>& constraints);

  template<typename Scalar, int XDim, int UDim>
  std::pair<std::vector<Vector<Scalar, XDim>>, std::vector<Vector<Scalar, UDim>>> getStateAndInputTrajectory(
    const std::vector<int>& numStates, const std::vector<int>& numInputs, const Vector<Scalar, Eigen::Dynamic>& w);

  /**
   * Extracts the problem state and inputs dimensions as well as number of constraints from a linear quadratic approximation
   * Looks at the size of the flowmap derivatives of the dynamics.
   * @return { numStatesPerStage, numInputsPerStage, numConstraintsPerStage}
   */
  template <typename Scalar, int XDim, int UDim>
  std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> getNumStatesInputsConstraints(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& linearQuadraticApproximation) {
    if (linearQuadraticApproximation.empty()) {
      return { std::vector<int>(0), std::vector<int>(0), std::vector<int>(0) };
    }

    const int N = linearQuadraticApproximation.size() - 1;
    std::vector<int> numStates;
    std::vector<int> numInputs;
    std::vector<int> numConstraints;
    numStates.reserve(N + 1);
    numInputs.reserve(N);
    numConstraints.reserve(N + 1);

    for (int k = 0; k < N; ++k) {
      numStates.push_back(linearQuadraticApproximation[k].dynamics.dfdx.cols());
      numInputs.push_back(linearQuadraticApproximation[k].dynamics.dfdu.cols());
      numConstraints.push_back(linearQuadraticApproximation[k].constraints.f.size());
    }
    numStates.push_back(linearQuadraticApproximation[N - 1].dynamics.dfdx.rows());
    numConstraints.push_back(linearQuadraticApproximation[N].constraints.f.size());

    return { numStates, numInputs, numConstraints };
  }

  /** Counts the number of decision variables in the QP */
  int getNumDecisionVariables(const std::vector<int>& numStates, const std::vector<int>& numInputs) {
    const auto totalNumberOfStates = std::accumulate(numStates.begin(), numStates.end(), 0);
    const auto totalNumberOfInputs = std::accumulate(numInputs.begin(), numInputs.end(), 0);
    return totalNumberOfStates + totalNumberOfInputs;
  }

  /** Counts the number of constraints in the QP */
  int getNumConstraints(const std::vector<int>& numStates, const std::vector<int>& numConstraints) {
    // Each stage constrains x_{k+1} states, adding the x_0 constraint, all states are constrained exactly once.
    return std::accumulate(numStates.begin(), numStates.end(), 0) + std::accumulate(numConstraints.begin(), numConstraints.end(), 0);
  }


  /**
   * Solves the discretized linear quadratic optimal control problem by constructing a dense QP and inverting the full KKT system.
   * The decision vector is defined as w = [dx[0], du[0], dx[1],  du[1], ..., dx[N]]
   *
   * @param lqApproximation : vector of stage-wise discrete quadratic cost and linear dynamics
   * @param dx0 : initial state deviation from the nominal trajectories.
   * @return trajectory of state and inputs (in relative coordinates), .i.e. dx(t), du(t)
   */
  template <typename Scalar, int XDim, int UDim>
  std::pair<std::vector<Vector<Scalar, XDim>>, std::vector<Vector<Scalar, UDim>>> solveLinearQuadraticProblem(const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqApproximation,
    const Vector<Scalar, XDim>& dx0)
  {
    // Extract sizes
    std::vector<int> numStates;
    std::vector<int> numInputs;
    std::vector<int> numConstraints;
    std::tie(numStates, numInputs, numConstraints) = getNumStatesInputsConstraints(lqApproximation);
    const auto numDecisionVariables = getNumDecisionVariables(numStates, numInputs);
    const auto numQpConstraints = getNumConstraints(numStates, numConstraints);

    // Construct QP
    const auto qpConstraints = getConstraintMatrices(lqApproximation, dx0, numQpConstraints, numDecisionVariables);
    const auto qpCosts = getCostMatrices(lqApproximation, numDecisionVariables);

    // Solve
    const auto primalDualSolution = solveDenseQp(qpCosts, qpConstraints);

    // Extract solution
    return getStateAndInputTrajectory<Scalar, XDim, UDim>(numStates, numInputs, primalDualSolution.first);
  }


  /**
   * Constructs the matrix of stacked dynamic constraints A w + b = 0
   *
   * A = [-I  *
   *       C  D  *  *
   *       A  B -I  *
   *       *  *  C  D  *  *
   *       *  *  A  B -I  *
   *       *  *  *  *  C  D  *
   *       *  *  *  *  A  B -I ]
   *       *  *  *  *  *  *  C ]
   *
   * b = [x0; e[0]; b[0]; ... e[N-1]; b[N-1]; e[N]]
   *
   * @param lqp : linear quadratic problem.
   * @param dx0 : initial state deviation from the nominal trajectories.
   * @param numConstraints : number of rows in A
   * @param numDecisionVariables : size of w
   * @return linear constraints in w, where w is the vector of decision variables
   */
  template <typename Scalar, int XDim, int UDim>
  VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0> getConstraintMatrices(const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp, const Vector<Scalar, XDim>& dx0,
    int numConstraints, int numDecisionVariables)
  {
    if (lqp.empty()) {
      return VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>();
    }

    const int N = lqp.size() - 1;

    // Preallocate full constraint matrix
    VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0> constraints;
    auto& A = constraints.dfdx;
    auto& b = constraints.f;
    A.setZero(numConstraints, numDecisionVariables);
    b.setZero(numConstraints);

    // Initial state constraint
    const int nx_0 = dx0.size();
    A.topLeftCorner(nx_0, nx_0) = -Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(nx_0, nx_0);
    b.topRows(nx_0) = dx0;

    int currRow = nx_0;
    int currCol = 0;
    for (int k = 0; k < N; ++k) {
      const auto& dynamics_k = lqp[k].dynamics;
      const auto& constraints_k = lqp[k].constraints;
      const int nu_k = dynamics_k.dfdu.cols();
      const int nx_k = dynamics_k.dfdx.cols();
      const int nx_Next = dynamics_k.dfdx.rows();
      const int nc_k = constraints_k.f.size();

      if (nc_k > 0) {
        // Add [C, D, 0]
        A.block(currRow, currCol, nc_k, nx_k + nu_k) << constraints_k.dfdx, constraints_k.dfdu;
        // Add [e]
        b.segment(currRow, nc_k) = constraints_k.f;

        currRow += nc_k;
      }

      // Add [A, B, -I]
      A.block(currRow, currCol, nx_Next, nx_k + nu_k + nx_Next) << dynamics_k.dfdx, dynamics_k.dfdu, -Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(nx_Next, nx_Next);
      // Add [b]
      b.segment(currRow, nx_Next) = dynamics_k.f;

      currRow += nx_Next;
      currCol += nx_k + nu_k;
    }

    // Final state constraint
    const auto& constraints_N = lqp[N].constraints;
    const int nc_N = constraints_N.f.size();
    if (nc_N > 0) {
      A.bottomRightCorner(nc_N, constraints_N.dfdx.cols()) = constraints_N.dfdx;
      b.bottomRows(nc_N) = constraints_N.f;
    }

    return constraints;
  }


  /**
   * Constructs a matrix of stacked cost functions  1/2 w' H w + g' w
   *
   * H = [ Q  P'
   *       P  R
   *             Q  P'
   *             P  R
   *                   Qf ]
   *
   * g = [q[0]; r[0]; q[1]; r[1]; ... qf ]
   *
   * @param lqp
   * @param numDecisionVariables : size of w
   * @return quadratic cost function in w, where w is the vector of decision variables
   */
  template <typename Scalar, int XDim, int UDim>
  ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0> getCostMatrices(const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp, int numDecisionVariables)
  {
    if (lqp.empty()) {
      return ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>();
    }

    const int N = lqp.size() - 1;

    // Preallocate full Cost matrices
    ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0> qpCost;
    auto& H = qpCost.dfdxx;
    auto& g = qpCost.dfdx;
    auto& c = qpCost.f;
    H.setZero(numDecisionVariables, numDecisionVariables);
    g.setZero(numDecisionVariables);
    c = 0.0;

    int currRow = 0;
    for (int k = 0; k < N; ++k) {
      const auto& cost_k = lqp[k].cost;
      const int nx_k = cost_k.dfdux.cols();
      const int nu_k = cost_k.dfdux.rows();

      // Add [ Q, P'
      //       P, R ]
      H.block(currRow, currRow, nx_k + nu_k, nx_k + nu_k) << cost_k.dfdxx, cost_k.dfdux.transpose(), cost_k.dfdux, cost_k.dfduu;
      // Add [ q, r]
      g.segment(currRow, nx_k + nu_k) << cost_k.dfdx, cost_k.dfdu;
      // Add nominal cost
      c += cost_k.f;

      currRow += nx_k + nu_k;
    }

    // Terminal cost
    const auto& cost_N = lqp[N].cost;
    const int nx_N = cost_N.dfdx.size();
    H.block(currRow, currRow, nx_N, nx_N) << cost_N.dfdxx;
    g.segment(currRow, nx_N) << cost_N.dfdx;
    c += cost_N.f;

    return qpCost;
  }

  /**
   * Constructs the equality constrained QP from the discretized linear quadratic control problem.
   * @param lqp : linear quadratic problem.
   * @param dx0 : initial state deviation from the nominal trajectories.
   * @return { cost matrices, constraint matrices }
   */
  template <typename Scalar, int XDim, int UDim>
  std::pair<ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>, VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>>
    getDenseQp(const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
      const Vector<Scalar, XDim>& dx0) {
    // Extract sizes
    std::vector<int> numStates;
    std::vector<int> numInputs;
    std::vector<int> numConstraints;
    std::tie(numStates, numInputs, numConstraints) = getNumStatesInputsConstraints(lqp);
    const auto numDecisionVariables = getNumDecisionVariables(numStates, numInputs);
    const auto numQpConstraints = getNumConstraints(numStates, numConstraints);

    // Construct QP
    const auto qpCosts = getCostMatrices(lqp, numDecisionVariables);
    const auto qpConstraints = getConstraintMatrices(lqp, dx0, numQpConstraints, numDecisionVariables);

    return { qpCosts, qpConstraints };
  }

  /**
   * Solves the equality constrained QP
   * min_w  1/2 w' H w + g' w
   *   s.t. A w + b = 0
   *
   *   Assumes H is positive definite, rows of A are linearly independent.
   *
   * @return {w, lambda} at the solution, where w is the vector of decision variables, and lambda is the vector of lagrange multipliers
   */
  template <typename Scalar>
  std::pair<Vector<Scalar, Eigen::Dynamic>, Vector<Scalar, Eigen::Dynamic>> solveDenseQp(const ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>& cost,
    const VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>& constraints) {
    const int m = constraints.dfdx.rows();
    const int n = constraints.dfdx.cols();

    // Assemble KKT condition
    Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> kktMatrix(n + m, n + m);
    Vector<Scalar, Eigen::Dynamic> kktRhs(n + m);
    kktMatrix << cost.dfdxx, constraints.dfdx.transpose(), constraints.dfdx,
      Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, m);
    kktRhs << -cost.dfdx, -constraints.f;

    // prerequisite for the LU factorization, and the solution would be non-unique.
    if (kktMatrix.fullPivLu().rank() != n + m) {
      throw std::runtime_error("KKT matrix is not full rank");
    }
    Vector<Scalar, Eigen::Dynamic>  sol = kktMatrix.lu().solve(kktRhs);
    return { sol.head(n), sol.tail(m) };
  }
  /**
   * Reconstructs the optimal state and input trajectory recursively based on the full qp solution vector
   * @param numStates : number of states per stage
   * @param numInputs : number of inputs per stage
   * @param w : the vector of decision variables
   * @return { state_trajectory, input_trajectory }
   */
  template<typename Scalar, int XDim, int UDim>
  std::pair<std::vector<Vector<Scalar, XDim>>, std::vector<Vector<Scalar, UDim>>>
    getStateAndInputTrajectory(const std::vector<int>& numStates, const std::vector<int>& numInputs,
      const Vector<Scalar, Eigen::Dynamic>& w)
  {
    assert(numStates.size() == numInputs.size() + 1);

    const int N = numInputs.size();

    std::vector<Vector<Scalar, XDim>> stateTrajectory;
    std::vector<Vector<Scalar, UDim>> inputTrajectory;
    stateTrajectory.reserve(N + 1);
    inputTrajectory.reserve(N);

    int index = 0;
    for (int k = 0; k < N; ++k) {
      // x[k]
      stateTrajectory.emplace_back(w.segment(index, XDim));
      index += XDim;

      // u[k]
      inputTrajectory.emplace_back(w.segment(index, UDim));
      index += UDim;
    }

    // x[N]
    stateTrajectory.emplace_back(w.segment(index, XDim));

    return { stateTrajectory, inputTrajectory };
  }
}  // namespace qp_solver
