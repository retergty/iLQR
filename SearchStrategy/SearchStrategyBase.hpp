#pragma once
#include "Types.hpp"
#include "SearchStrategySettings.hpp"
#include "Controller.hpp"
#include "LinearController.hpp"
#include "DualSolution.hpp"
#include "PrimalSolution.hpp"
#include "ProblemMetrics.hpp"
#include "PerformanceIndex.hpp"
#include "ModelData.hpp"

template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
struct SearchStrategySolution
{
  using PrimalSolution_t = PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength>;
  using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
  using ProblemMetrics_t = ProblemMetrics<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;

  Scalar avgTimeStep;
  DualSolution_t dualSolution;
  PrimalSolution_t primalSolution;
  ProblemMetrics_t problemMetrics;
  PerformanceIndex_t performanceIndex;
};

template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
struct SearchStrategySolutionRef
{
  using PrimalSolution_t = PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength>;
  using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
  using ProblemMetrics_t = ProblemMetrics<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using SearchStrategySolution_t = SearchStrategySolution<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;

  SearchStrategySolutionRef(SearchStrategySolution_t &s)
      : avgTimeStep(s.avgTimeStep),
        dualSolution(s.dualSolution),
        primalSolution(s.primalSolution),
        problemMetrics(s.problemMetrics),
        performanceIndex(s.performanceIndex)
  {
  }

  SearchStrategySolutionRef(Scalar avgTimeStepArg, DualSolution_t &dualSolutionArg, PrimalSolution_t &primalSolutionArg, ProblemMetrics_t &problemMetricsArg,
                            PerformanceIndex_t &performanceIndexArg)
      : avgTimeStep(avgTimeStepArg),
        dualSolution(dualSolutionArg),
        primalSolution(primalSolutionArg),
        problemMetrics(problemMetricsArg),
        performanceIndex(performanceIndexArg)
  {
  }

  Scalar &avgTimeStep;
  DualSolution_t &dualSolution;
  PrimalSolution_t &primalSolution;
  ProblemMetrics_t &problemMetrics;
  PerformanceIndex_t &performanceIndex;
};

/**
 * This class is an interface class for search strategies such as line-search, trust-region.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains, int StateIneqConstrains , int StateInputEqConstrains, int StateInputIneqConstrains,
          int FinalStateEqConstrains, int FinalStateIneqConstrains>
class SearchStrategyBase
{
public:
  using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
  using LinearController_t = LinearController<Scalar, XDimisions, UDimisions, PredictLength + 1>;
  using StateVector_t = Vector<Scalar, XDimisions>;
  using PerformanceIndex_t = PerformanceIndex<Scalar>;
  using ModelData_t = ModelData<Scalar, XDimisions, UDimisions>;
  using SearchStrategySolution_t = SearchStrategySolution<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
  using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;

  /**
   * Constructor.
   * @param [in] baseSettings: The basic settings for the search strategy algorithms.
   */
  explicit SearchStrategyBase() {}

  virtual ~SearchStrategyBase() = default;
  SearchStrategyBase(const SearchStrategyBase &) = delete;
  SearchStrategyBase &operator=(const SearchStrategyBase &) = delete;

  /**
   * Resets the class to its state after construction.
   */
  virtual void reset() = 0;

  /**
   * Finds the optimal trajectories, controller, and performance index based on the given controller and its increment.
   *
   * @param [in] initialTime: Initial time.
   * @param [in] finalTime: final time
   * @param [in] initState: Initial state
   * @param [in] expectedCost: The expected cost based on the LQ model optimization.
   * @param [in] unoptimizedController: The unoptimized controller which search will be performed.
   * @param [in] dualSolution: The dual solution.
   * @param [in/out]
   * @param [out] solution: Output of search (primalSolution, performanceIndex, problemMetrics, avgTimeStep)
   * @return whether the search was successful or failed.
   */
  virtual bool run(const std::pair<Scalar, Scalar> &timePeriod, const StateVector_t &initState, const Scalar expectedCost,
                   const LinearController_t &unoptimizedController, const DualSolution_t &dualSolution, SearchStrategySolutionRef_t &solution) = 0;

  /**
   * Checks convergence of the main loop of DDP.
   *
   * @param [in] unreliableControllerIncrement: True if the controller is designed based on an unreliable LQ approximation such as
   * operating trajectories
   * @param [in] previousPerformanceIndex: The previous iteration's PerformanceIndex.
   * @param [in] currentPerformanceIndex: The current iteration's PerformanceIndex.
   * @return A pair of (isOptimizationConverged, infoString)
   */
  virtual bool checkConvergence(bool unreliableControllerIncrement,
                                const PerformanceIndex_t &previousPerformanceIndex,
                                const PerformanceIndex_t &currentPerformanceIndex) const = 0;

  /**
   * Computes the Riccati modification based on the strategy.
   *
   * @param [in] projectedModelData: The projected data model
   * @param [out] deltaQm: The Riccati modifier to cost 2nd derivative w.r.t. state.
   * @param [out] deltaGv: The Riccati modifier to cost derivative w.r.t. input.
   * @param [out] deltaGm: The Riccati modifier to cost input-state derivative.
   */
  virtual void computeRiccatiModification(const ModelData_t &projectedModelData, Matrix<Scalar, XDimisions, XDimisions> &deltaQm) const = 0;

  /**
   * Augments the Hessian of Hamiltonian based on the strategy.
   *
   * @param [in] modelData: The model data.
   * @param [in] Hm: The Hessian of Hamiltonian that should be augmented.
   * @return The augmented Hamiltonian's Hessian.
   */
  virtual Matrix<Scalar, UDimisions, UDimisions> augmentHamiltonianHessian(const ModelData_t &modelData, const Matrix<Scalar, UDimisions, UDimisions> &Hm) const = 0;

protected:
  constexpr static SearchStrategyBaseSettings<Scalar> baseSettings_{};
};