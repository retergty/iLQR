/**
 * @file SearchStrategyBase.hpp
 * @brief 搜索策略基类：线搜索/置信域等子问题求解接口及搜索解容器。
 */
#pragma once
#include <utility>

template <typename Descriptor>
struct iLQRTypes;
/**
 * @brief 搜索策略外部候选缓冲区视图：绑定 baseline/candidate 两个候选槽。
 *
 * 该结构不拥有 primalSolution/problemMetrics 存储，只把线搜索需要反复覆盖的
 * 两组外部缓冲区组织成固定双槽。搜索结束后，selectedSlotPtr 指向被接受的槽。
 *
 * @tparam Descriptor iLQR 类型描述。
 */
template <typename Descriptor>
struct SearchStrategySolutionBuffer {
  using Types = iLQRTypes<Descriptor>;
  using PrimalSolution_t = typename Types::PrimalSolution_t;
  using ProblemMetrics_t = typename Types::ProblemMetrics_t;
  using PerformanceIndex_t = typename Types::PerformanceIndex_t;

  /** @brief 单个候选槽：轨迹/metrics 绑定外部存储，performanceIndex 本地保存。
   */
  struct Slot {
    PrimalSolution_t& primalSolution;
    ProblemMetrics_t& problemMetrics;
    PerformanceIndex_t performanceIndex;
  };

  SearchStrategySolutionBuffer(PrimalSolution_t& baselinePrimalSolution,
                               ProblemMetrics_t& baselineProblemMetrics,
                               PrimalSolution_t& candidatePrimalSolution,
                               ProblemMetrics_t& candidateProblemMetrics)
      : baseline{baselinePrimalSolution, baselineProblemMetrics, {}},
        candidate{candidatePrimalSolution, candidateProblemMetrics, {}},
        selectedSlotPtr(&baseline) {}

  /** @brief 零步长候选，作为 Armijo 比较基准。 */
  Slot baseline;
  /** @brief 正步长候选，线搜索过程中会被重复覆盖。 */
  Slot candidate;
  /** @brief 搜索结束后指向被接受的候选槽。 */
  Slot* selectedSlotPtr;
};

/**
 * @brief 搜索策略抽象基类：线搜索、置信域等子问题求解的通用接口。
 * @tparam Descriptor iLQR 类型描述。
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
  using SearchStrategySolutionBuffer_t =
      SearchStrategySolutionBuffer<Descriptor>;
  using SearchStrategySolutionSlot_t =
      typename SearchStrategySolutionBuffer_t::Slot;
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
   * @brief 根据未优化控制器与对偶解执行搜索，在候选缓冲区中返回被选中的槽。
   * @param [in] timePeriod 时间区间 (初始时间, 终止时间)。
   * @param [in] initState 初始状态。
   * @param [in] expectedCost 基于 LQ 模型的期望代价。
   * @param [in,out] solutionBuffer 外部 baseline/candidate 候选缓冲区。
   * @param [in] unoptimizedController 待搜索的未优化控制器。
   * @param [in] dualSolution 对偶解。
   * @return 被选中的候选槽；返回 nullptr 表示搜索失败且没有可用候选。
   */
  virtual SearchStrategySolutionSlot_t* run(
      const std::pair<Scalar, Scalar>& timePeriod,
      const StateVector_t& initState, const Scalar expectedCost,
      SearchStrategySolutionBuffer_t& solutionBuffer,
      const LinearController_t& unoptimizedController,
      const DualSolution_t& dualSolution) = 0;

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