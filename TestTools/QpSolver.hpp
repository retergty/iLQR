#pragma once

#include <Eigen/LU>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "QpSolverTypes.hpp"
#include "QpTrajectories.hpp"

namespace qp_solver {

template <typename Scalar, int XDim, int UDim>
VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>
getConstraintMatrices(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
    const Vector<Scalar, XDim>& dx0, int numConstraints,
    int numDecisionVariables);

template <typename Scalar, int XDim, int UDim>
ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0> getCostMatrices(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
    int numDecisionVariables);

template <typename Scalar>
std::pair<Vector<Scalar, Eigen::Dynamic>, Vector<Scalar, Eigen::Dynamic>>
solveDenseQp(
    const ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>& cost,
    const VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic,
                                            Eigen::Dynamic, 0>& constraints);

template <typename Scalar, int XDim, int UDim>
std::pair<std::vector<Vector<Scalar, XDim>>, std::vector<Vector<Scalar, UDim>>>
getStateAndInputTrajectory(const std::vector<int>& numStates,
                           const std::vector<int>& numInputs,
                           const Vector<Scalar, Eigen::Dynamic>& w);

/**
 * 提取问题的状态、输入维度以及约束数量。
 * 这些信息来自线性二次近似，并根据动力学 flowmap 导数尺寸确定。
 * 动力学 flowmap 导数尺寸。
 * @return { numStatesPerStage, numInputsPerStage, numConstraintsPerStage }。
 */
template <typename Scalar, int XDim, int UDim>
std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
getNumStatesInputsConstraints(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>&
        linearQuadraticApproximation) {
  if (linearQuadraticApproximation.empty()) {
    return {std::vector<int>(0), std::vector<int>(0), std::vector<int>(0)};
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
    numConstraints.push_back(
        linearQuadraticApproximation[k].constraints.f.size());
  }
  numStates.push_back(linearQuadraticApproximation[N - 1].dynamics.dfdx.rows());
  numConstraints.push_back(
      linearQuadraticApproximation[N].constraints.f.size());

  return {numStates, numInputs, numConstraints};
}

/** 统计 QP 中的决策变量数量。 */
int getNumDecisionVariables(const std::vector<int>& numStates,
                            const std::vector<int>& numInputs) {
  const auto totalNumberOfStates =
      std::accumulate(numStates.begin(), numStates.end(), 0);
  const auto totalNumberOfInputs =
      std::accumulate(numInputs.begin(), numInputs.end(), 0);
  return totalNumberOfStates + totalNumberOfInputs;
}

/** 统计 QP 中的约束数量。 */
int getNumConstraints(const std::vector<int>& numStates,
                      const std::vector<int>& numConstraints) {
  // 每个阶段约束 x_{k+1} 状态；加上 x_0 约束后，所有状态
  // 都会被精确约束一次。
  return std::accumulate(numStates.begin(), numStates.end(), 0) +
         std::accumulate(numConstraints.begin(), numConstraints.end(), 0);
}

/**
 * 通过构造稠密 QP 并反转完整 KKT 系统，求解离散化线性二次最优控制问题。
 * 构造稠密 QP 并反转完整 KKT 系统。决策
 * 向量定义为 w = [dx[0], du[0], dx[1], du[1], ..., dx[N]]。
 *
 * @param lqApproximation: 按阶段存储的离散二次代价和线性动力学向量。
 * @param dx0: 相对于名义轨迹的初始状态偏差。
 * @return 状态和输入轨迹（相对坐标），即。
 * dx(t), du(t)
 */
template <typename Scalar, int XDim, int UDim>
std::pair<std::vector<Vector<Scalar, XDim>>, std::vector<Vector<Scalar, UDim>>>
solveLinearQuadraticProblem(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>&
        lqApproximation,
    const Vector<Scalar, XDim>& dx0) {
  // 提取尺寸
  std::vector<int> numStates;
  std::vector<int> numInputs;
  std::vector<int> numConstraints;
  std::tie(numStates, numInputs, numConstraints) =
      getNumStatesInputsConstraints(lqApproximation);
  const auto numDecisionVariables =
      getNumDecisionVariables(numStates, numInputs);
  const auto numQpConstraints = getNumConstraints(numStates, numConstraints);

  // 构造 QP
  const auto qpConstraints = getConstraintMatrices(
      lqApproximation, dx0, numQpConstraints, numDecisionVariables);
  const auto qpCosts = getCostMatrices(lqApproximation, numDecisionVariables);

  // 求解
  const auto primalDualSolution = solveDenseQp(qpCosts, qpConstraints);

  // 提取解。
  return getStateAndInputTrajectory<Scalar, XDim, UDim>(
      numStates, numInputs, primalDualSolution.first);
}

/**
 * 构造堆叠动力学约束矩阵 A w + b = 0。
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
 * @param lqp: 线性二次问题。
 * @param dx0: 相对于名义轨迹的初始状态偏差。
 * @param numConstraints: A 的行数。
 * @param numDecisionVariables: w 的尺寸。
 * @return w 上的线性约束，其中 w 是决策变量向量。
 */
template <typename Scalar, int XDim, int UDim>
VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>
getConstraintMatrices(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
    const Vector<Scalar, XDim>& dx0, int numConstraints,
    int numDecisionVariables) {
  if (lqp.empty()) {
    return VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic,
                                             Eigen::Dynamic, 0>();
  }

  const int N = lqp.size() - 1;

  // 预分配完整约束矩阵。
  VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic, Eigen::Dynamic, 0>
      constraints;
  auto& A = constraints.dfdx;
  auto& b = constraints.f;
  A.setZero(numConstraints, numDecisionVariables);
  b.setZero(numConstraints);

  // 初始状态约束。
  const int nx_0 = dx0.size();
  A.topLeftCorner(nx_0, nx_0) =
      -Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(nx_0, nx_0);
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
      // 添加 [C, D, 0]。
      A.block(currRow, currCol, nc_k, nx_k + nu_k) << constraints_k.dfdx,
          constraints_k.dfdu;
      // 添加 [e]。
      b.segment(currRow, nc_k) = constraints_k.f;

      currRow += nc_k;
    }

    // 添加 [A, B, -I]。
    A.block(currRow, currCol, nx_Next, nx_k + nu_k + nx_Next)
        << dynamics_k.dfdx,
        dynamics_k.dfdu,
        -Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(nx_Next,
                                                                  nx_Next);
    // 添加 [b]。
    b.segment(currRow, nx_Next) = dynamics_k.f;

    currRow += nx_Next;
    currCol += nx_k + nu_k;
  }

  // 终端状态约束。
  const auto& constraints_N = lqp[N].constraints;
  const int nc_N = constraints_N.f.size();
  if (nc_N > 0) {
    A.bottomRightCorner(nc_N, constraints_N.dfdx.cols()) = constraints_N.dfdx;
    b.bottomRows(nc_N) = constraints_N.f;
  }

  return constraints;
}

/**
 * 构造堆叠代价函数矩阵 1/2 w' H w + g' w。
 *
 * H = [ Q  P'
 *       P  R
 *             Q  P'
 *             P  R
 *                   Qf ]
 *
 * g = [q[0]; r[0]; q[1]; r[1]; ... qf ]
 *
 * @param lqp 线性二次问题。
 * @param numDecisionVariables: w 的尺寸。
 * @return w 上的二次代价函数，其中 w 是决策
 * 变量。
 */
template <typename Scalar, int XDim, int UDim>
ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0> getCostMatrices(
    const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
    int numDecisionVariables) {
  if (lqp.empty()) {
    return ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>();
  }

  const int N = lqp.size() - 1;

  // 预分配完整代价矩阵。
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

    // 添加 [ Q, P'。
    //       P, R ]
    H.block(currRow, currRow, nx_k + nu_k, nx_k + nu_k) << cost_k.dfdxx,
        cost_k.dfdux.transpose(), cost_k.dfdux, cost_k.dfduu;
    // 添加 [ q, r]。
    g.segment(currRow, nx_k + nu_k) << cost_k.dfdx, cost_k.dfdu;
    // 添加名义代价。
    c += cost_k.f;

    currRow += nx_k + nu_k;
  }

  // 终端代价。
  const auto& cost_N = lqp[N].cost;
  const int nx_N = cost_N.dfdx.size();
  H.block(currRow, currRow, nx_N, nx_N) << cost_N.dfdxx;
  g.segment(currRow, nx_N) << cost_N.dfdx;
  c += cost_N.f;

  return qpCost;
}

/**
 * 由离散化线性二次控制问题构造等式约束 QP。
 * 控制问题。
 * @param lqp: 线性二次问题。
 * @param dx0: 相对于名义轨迹的初始状态偏差。
 * @return { 代价矩阵, 约束矩阵 }。
 */
template <typename Scalar, int XDim, int UDim>
std::pair<ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>,
          VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic,
                                            Eigen::Dynamic, 0>>
getDenseQp(const std::vector<LinearQuadraticStage<Scalar, XDim, UDim>>& lqp,
           const Vector<Scalar, XDim>& dx0) {
  // 提取尺寸
  std::vector<int> numStates;
  std::vector<int> numInputs;
  std::vector<int> numConstraints;
  std::tie(numStates, numInputs, numConstraints) =
      getNumStatesInputsConstraints(lqp);
  const auto numDecisionVariables =
      getNumDecisionVariables(numStates, numInputs);
  const auto numQpConstraints = getNumConstraints(numStates, numConstraints);

  // 构造 QP
  const auto qpCosts = getCostMatrices(lqp, numDecisionVariables);
  const auto qpConstraints =
      getConstraintMatrices(lqp, dx0, numQpConstraints, numDecisionVariables);

  return {qpCosts, qpConstraints};
}

/**
 * 求解等式约束 QP。
 * min_w  1/2 w' H w + g' w。
 *   s.t. A w + b = 0
 *
 *   假定 H 正定且 A 的各行线性无关。
 *
 * @return 解处的 {w, lambda}，其中 w 是决策
 * 变量向量，lambda 是拉格朗日乘子向量。
 */
template <typename Scalar>
std::pair<Vector<Scalar, Eigen::Dynamic>, Vector<Scalar, Eigen::Dynamic>>
solveDenseQp(
    const ScalarFunctionQuadraticApproximation<Scalar, Eigen::Dynamic, 0>& cost,
    const VectorFunctionLinearApproximation<Scalar, Eigen::Dynamic,
                                            Eigen::Dynamic, 0>& constraints) {
  const int m = constraints.dfdx.rows();
  const int n = constraints.dfdx.cols();

  // 组装 KKT 条件。
  Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> kktMatrix(n + m, n + m);
  Vector<Scalar, Eigen::Dynamic> kktRhs(n + m);
  kktMatrix << cost.dfdxx, constraints.dfdx.transpose(), constraints.dfdx,
      Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, m);
  kktRhs << -cost.dfdx, -constraints.f;

  // LU 分解的前提条件，否则解会
  // 不唯一。
  if (kktMatrix.fullPivLu().rank() != n + m) {
    throw std::runtime_error("KKT matrix is not full rank");
  }
  Vector<Scalar, Eigen::Dynamic> sol = kktMatrix.lu().solve(kktRhs);
  return {sol.head(n), sol.tail(m)};
}
/**
 * 基于完整 QP 解向量递归重建最优状态和输入轨迹。
 * 完整 QP 解向量。
 * @param numStates: 每阶段状态数量。
 * @param numInputs: 每阶段输入数量。
 * @param w: 决策变量向量。
 * @return { state_trajectory, input_trajectory }。
 */
template <typename Scalar, int XDim, int UDim>
std::pair<std::vector<Vector<Scalar, XDim>>, std::vector<Vector<Scalar, UDim>>>
getStateAndInputTrajectory(const std::vector<int>& numStates,
                           const std::vector<int>& numInputs,
                           const Vector<Scalar, Eigen::Dynamic>& w) {
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

  return {stateTrajectory, inputTrajectory};
}
}  // namespace qp_solver
