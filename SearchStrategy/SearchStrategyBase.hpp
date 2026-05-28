/**
 * @file SearchStrategyBase.hpp
 * @brief 搜索策略基类：线搜索/置信域等子问题求解接口及搜索解容器。
 */
#pragma once
#include "Matrix/Types.hpp"

template <typename Descriptor>
struct iLQRTypes;
/**
 * @brief 搜索策略候选输出：平均步长、原始解、问题指标与性能指标。
 * @note 候选解不缓存对偶解；对偶解由 SearchStrategySolutionRef 在外部保存。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam PredictLength 预测步数。
 * @tparam StateEqConstrains 等 约束维度。
 */
template <typename Descriptor>
struct SearchStrategySolution {
  using Types = iLQRTypes<Descriptor>;
  using Scalar = typename Types::Scalar;
  using PrimalSolution_t = typename Types::PrimalSolution_t;
  using DualSolution_t = typename Types::DualSolution_t;
  using ProblemMetrics_t = typename Types::ProblemMetrics_t;
  using PerformanceIndex_t = typename Types::PerformanceIndex_t;

  /** @brief 平均时间步长。 */
  Scalar avgTimeStep;
  /** @brief 原始解。 */
  PrimalSolution_t primalSolution;
  /** @brief 问题指标。 */
  ProblemMetrics_t problemMetrics;
  /** @brief 性能指标。 */
  PerformanceIndex_t performanceIndex;
};

/** @brief 搜索策略写回视图：绑定外部 avg/dual/primal/metrics/performance。 */
template <typename Descriptor>
struct SearchStrategySolutionRef {
  using Types = iLQRTypes<Descriptor>;
  using Scalar = typename Types::Scalar;
  using PrimalSolution_t = typename Types::PrimalSolution_t;
  using DualSolution_t = typename Types::DualSolution_t;
  using ProblemMetrics_t = typename Types::ProblemMetrics_t;
  using PerformanceIndex_t = typename Types::PerformanceIndex_t;

  using SearchStrategySolution_t = SearchStrategySolution<Descriptor>;

  /** @brief 直接绑定各成员引用。 */
  SearchStrategySolutionRef(Scalar& avgTimeStepArg,
                            DualSolution_t& dualSolutionArg,
                            PrimalSolution_t& primalSolutionArg,
                            ProblemMetrics_t& problemMetricsArg,
                            PerformanceIndex_t& performanceIndexArg)
      : avgTimeStep(avgTimeStepArg),
        dualSolution(dualSolutionArg),
        primalSolution(primalSolutionArg),
        problemMetrics(problemMetricsArg),
        performanceIndex(performanceIndexArg) {}

  /** @brief 平均时间步长引用。 */
  Scalar& avgTimeStep;
  /** @brief 对偶解引用。 */
  DualSolution_t& dualSolution;
  /** @brief 原始解引用。 */
  PrimalSolution_t& primalSolution;
  /** @brief 问题指标引用。 */
  ProblemMetrics_t& problemMetrics;
  /** @brief 性能指标引用。 */
  PerformanceIndex_t& performanceIndex;
};

/**
 * @brief 搜索策略抽象基类：线搜索、置信域等子问题求解的通用接口。
 * @tparam Scalar 标量类型。
 * @tparam XDim 状态维度。
 * @tparam UDim 控制维度。
 * @tparam PredictLength 预测步数。
 * @tparam StateEqConstrains 等 约束维度。
 */
template <typename Descriptor>
class SearchStrategyBase {
 public:
  using Types = iLQRTypes<Descriptor>;
  using Scalar = typename Types::Scalar;
  using LinearController_t = typename Types::LinearController_t;
  using StateVector_t = typename Types::StateVector_t;
  using PerformanceIndex_t = typename Types::PerformanceIndex_t;
  using ModelData_t = typename Types::ModelData_t;
  using DualSolution_t = typename Types::DualSolution_t;
  using SearchStrategySolution_t = SearchStrategySolution<Descriptor>;
  using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Descriptor>;
  using HmMatrix_t = typename Types::HmMatrix_t;
  using SmMatrix_t = typename Types::SmMatrix_t;

  /** @brief 默认构造。 */
  explicit SearchStrategyBase() {}

  /** @brief 虚析构。 */
  virtual ~SearchStrategyBase() = default;
  SearchStrategyBase(const SearchStrategyBase&) = delete;
  SearchStrategyBase& operator=(const SearchStrategyBase&) = delete;

  /** @brief 重置为构造后状态。 */
  virtual void reset() = 0;

  /**
   * @brief 根据未优化控制器与对偶解执行搜索，将最优轨迹、控制器与性能指标写入
   * 解。
   * @param [in] timePeriod 时间区间 (初始时间, 终止时间)。
   * @param [in] initState 初始状态。
   * @param [in] expectedCost 基于 LQ 模型的期望代价。
   * @param [in] unoptimizedController 待搜索的未优化控制器。
   * @param [in] dualSolution 对偶解。
   * @param [in,out] solution 解。
   * 输出（primalSolution、performanceIndex、problemMetrics、avgTimeStep）。
   * @return 搜索是否成功。
   */
  virtual bool run(const std::pair<Scalar, Scalar>& timePeriod,
                   const StateVector_t& initState, const Scalar expectedCost,
                   const LinearController_t& unoptimizedController,
                   const DualSolution_t& dualSolution,
                   SearchStrategySolutionRef_t& solution) = 0;

  /**
   * @brief 检查 DDP 主循环是否收敛。
   * @param [in] unreliableControllerIncrement 控制器是否基于不可靠 LQ
   * 近似（如工作点轨迹）设计。
   * @param [in] previousPerformanceIndex 上一迭代性能指标。
   * @param [in] currentPerformanceIndex 当前迭代性能指标。
   * @return 是否已收敛。
   */
  virtual bool checkConvergence(
      bool unreliableControllerIncrement,
      const PerformanceIndex_t& previousPerformanceIndex,
      const PerformanceIndex_t& currentPerformanceIndex) const = 0;

  /**
   * @brief 根据策略计算 Riccati 修正（如代价对状态的二阶修正 deltaQm）。
   * @param [in] modelData 当前节点模型数据。
   * @param [out] deltaQm 代价对状态二阶导的 Riccati 修正。
   */
  virtual void computeRiccatiModification(const ModelData_t& modelData,
                                          SmMatrix_t& deltaQm) const = 0;

  /**
   * @brief 根据策略对哈密顿量 Hessian 进行增广（如数值稳定性修正）。
   * @param [in] modelData 模型数据。
   * @param [in] Hm 待增广的哈密顿量 Hessian。
   * @return 增广后的哈密顿量 Hessian。
   */
  virtual HmMatrix_t augmentHamiltonianHessian(const ModelData_t& modelData,
                                               const HmMatrix_t& Hm) const = 0;
};