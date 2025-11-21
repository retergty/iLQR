#pragma once
#include "Types.hpp"
#include "Dynamics.hpp"
#include "Cost.hpp"
#include "Controller.hpp"
#include "OptimalControlProblem.hpp"
#include "SensitivityIntegrator.hpp"
#include "PrimalSolution.hpp"
#include "DualSolution.hpp"
#include "LinearQuadraticApproximator.hpp"
#include "DDPData.hpp"
#include "ProblemMetrics.hpp"
#include "RolloutBase.hpp"
#include "DDPHelperFunction.hpp"
#include "PerformanceIndex.hpp"
#include "HessianCorrection.hpp"
#include "DDPSetting.hpp"
#include "LinearController.hpp"
#include "SearchStrategyBase.hpp"
#include "LineSearchStrategy.hpp"
#include "DiscreteTimeRiccatiEquations.hpp"
#include "LinearAlgebra.hpp"
#include "OptimalControlProblemHelperFunction.hpp"
#include "ChangeOfInputVariables.hpp"
#include "TimeTriggeredRollout.hpp"

template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains = 0, int StateIneqConstrains = 0, int StateInputEqConstrains = 0, int StateInputIneqConstrains = 0,
          int FinalStateEqConstrains = 0, int FinalStateIneqConstrains = 0>
class iLQR
{
public:
    using OptimalControlProblem_t = OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using ControlledSystemBase_t = ControlledSystemBase<Scalar, XDimisions, UDimisions>;
    using StateVector_t = Vector<Scalar, XDimisions>;
    using InputVector_t = Vector<Scalar, UDimisions>;
    using LvVector_t = Vector<Scalar, UDimisions>;
    using KmMatrix_t = Matrix<Scalar, UDimisions, XDimisions>;
    using SmMatrix_t = Matrix<Scalar, XDimisions, XDimisions>;
    using SvVector_t = Vector<Scalar, XDimisions>;
    using GmMatrix_t = Matrix<Scalar, UDimisions, XDimisions>;
    using HmMatrix_t = Matrix<Scalar, UDimisions, UDimisions>;
    using GvVector_t = Vector<Scalar, UDimisions>;

    using ModelData_t = ModelData<Scalar, XDimisions, UDimisions>;
    using IntermediateMultiplierCollection_t = MultiplierCollection<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>;
    using FinalMultiplierCollection_t = MultiplierCollection<Scalar, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0>;
    using IntermediateMetrics_t = Metrics<Scalar, XDimisions, UDimisions, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains>;
    using FinalMetrics_t = Metrics<Scalar, XDimisions, UDimisions, FinalStateEqConstrains, FinalStateIneqConstrains, 0, 0>;
    using RolloutBase_t = RolloutBase<Scalar, XDimisions, UDimisions, PredictLength>;
    using TimeTriggeredRollout_t = TimeTriggeredRollout<Scalar, XDimisions, UDimisions, PredictLength>;
    using SearchStrategySolution_t = SearchStrategySolution<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;

    using LinearController_t = LinearController<Scalar, XDimisions, UDimisions, PredictLength>;
    using SearchStrategyBase_t = SearchStrategyBase<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using RiccatiModification_t = RiccatiModification<Scalar, XDimisions, UDimisions>;

    using TimeTrajectory_t = std::array<Scalar, PredictLength>;
    using StateTrajectory_t = std::array<Vector<Scalar, XDimisions>, PredictLength>;
    using InputTrajectory_t = std::array<Vector<Scalar, UDimisions>, PredictLength>;
    using IntermediateMultiplierTrajectory_t = std::array<IntermediateMultiplierCollection_t, PredictLength>;
    using ModelDataTrajectory_t = std::array<ModelData_t, PredictLength>;

    using PrimalSolution_t = PrimalSolution<Scalar, XDimisions, UDimisions, PredictLength>;
    using DualSolution_t = DualSolution<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
    using DualSolutionRef_t = DualSolutionRef<Scalar, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains, PredictLength>;
    using LinearQuadraticApproximator_t = LinearQuadraticApproximator<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using PrimalDataContainer_t = PrimalDataContainer<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using DualDataContainer_t = DualDataContainer<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using ProblemMetrics_t = ProblemMetrics<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using PerformanceIndex_t = PerformanceIndex<Scalar>;

    using EK2DynamicsDiscretizer_t = EK2DynamicsDiscretizer<Scalar, XDimisions, UDimisions>;

    using ValueFunctionQuadraticApproximation_t = ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>;
    using DiscreteTimeRiccatiEquations_t = DiscreteTimeRiccatiEquations<Scalar, XDimisions, UDimisions>;

    iLQR(ControlledSystemBase_t *systemPtr, const Scalar timestep) : rollout_(systemPtr, timestep), predictTimeStep_(timestep) {

                                                                     };

    /**
     * The main routine of solver which runs the optimizer for a given initial state, initial time, and final time.
     *
     * @param [in] initTime: The initial time.
     * @param [in] initState: The initial state.
     */
    void run(Scalar initTime, const StateVector_t &initState)
    {
        // initialize parameters
        initTime_ = initTime;
        initState_ = initState;
        finalTime_ = initTime + predictTimeStep_ * (PredictLength);
        const auto initIteration = totalNumIterations_;

        // optimized --> nominal: initializes the nominal primal and dual solutions based on the optimized ones
        bool initialSolutionExists = initializePrimalSolution(); // true if the rollout is not purely from the Initializer
        initializeDualSolutionAndMetrics();

        // convergence variables of the main loop
        bool isConverged = false;

        // DDP main loop
        while (true)
        {
            // nominal --> nominal: constructs the LQ problem around the nominal trajectories
            approximateOptimalControlProblem();

            // nominal --> nominal: solves the LQ problem
            avgTimeStepBP_ = solveSequentialRiccatiEquations(nominalPrimalData_.modelDataFinalTime.cost);

            // calculate controller and store the result in unoptimizedController_
            calculateController();

            // the expected cost/merit calculated by the Riccati solution is not reliable
            const auto lqModelExpectedCost = initialSolutionExists ? nominalDualData_.valueFunctionTrajectory.front().f : performanceIndex_.merit;

            // nominal --> optimized: based on the current LQ solution updates the optimized primal and dual solutions
            takePrimalDualStep(lqModelExpectedCost);

            // iteration info
            ++totalNumIterations_;

            // check convergence
            isConverged = searchStrategyPtr_->checkConvergence(
                !initialSolutionExists, performanceIndexLast_, performanceIndex_);
            initialSolutionExists = true;

            if (isConverged || (totalNumIterations_ - initIteration) == ddpSettings_.maxNumIterations_)
            {
                break;
            }
            else
            {
                // optimized --> nominal: use the optimized solution as the nominal for the next iteration
                optimizedDualSolution_.swap(nominalDualData_.dualSolution);
                optimizedPrimalSolution_.swap(nominalPrimalData_.primalSolution);
                optimizedProblemMetrics_.swap(nominalPrimalData_.problemMetrics);
            }
        } // end of while loop
    }

private:
    /** Initializes the nominal primal based on the optimized ones.
     * @return True if the rollout is not purely from the Initializer.
     */
    bool initializePrimalSolution()
    {
        // try to initialize with controller
        bool initialSolutionExists = rolloutInitialController(optimizedPrimalSolution_, nominalPrimalData_.primalSolution);

        // if rolloutInitialController failed, try to initialize with PrimalSolution's state-input trajectories
        if (!initialSolutionExists)
        {
            // // display
            // if (ddpSettings_.displayInfo_) {
            //     std::cerr << "Initial controller is unavailable. Solver resorts to use PrimalSolution trajectories ...\n";
            // }

            initialSolutionExists = extractInitialTrajectories(optimizedPrimalSolution_, nominalPrimalData_.primalSolution);
        }
        assert(initialSolutionExists);

        return initialSolutionExists;
    }

    /**
     * Extracts the PrimalSolution trajectories from inputPrimalSolution. In general, it will try to extract in time period
     * [initTime, finalTime]. However, if inputPrimalSolution's timeTrajectory does not cover the period [initTime, finalTime],
     * it will extract the solution until the last time of the timeTrajectory
     *
     * @param [in] inputPrimalSolution: Its controller will be used for rollout.
     * @param [out] outputPrimalSolution: The resulting PrimalSolution.
     * @return True if the extraction was successful.
     */
    bool extractInitialTrajectories(PrimalSolution_t &inputPrimalSolution, PrimalSolution_t &outputPrimalSolution)
    {
        if (inputPrimalSolution.timeTrajectory_.empty())
        {
            return false;
        }

        // const auto finalTime = std::max(initTime_, std::min(inputPrimalSolution.timeTrajectory_.back(), finalTime_));

        if (initTime_ < finalTime_)
        {
            // if (ddpSettings_.debugPrintRollout_) {
            //     std::cerr << "[GaussNewtonDDP::extractInitialTrajectories] for t = [" << initTime_ << ", " << finalTime_ << "]\n";
            //     std::cerr << "\twill use PrimalSolution trajectory for t = [" << initTime_ << ", " << finalTime << "]\n";
            // }
            // extractPrimalSolution({ initTime_, finalTime_ }, inputPrimalSolution, outputPrimalSolution);

            // copy
            outputPrimalSolution.timeTrajectory_ = inputPrimalSolution.timeTrajectory_;
            outputPrimalSolution.stateTrajectory_ = inputPrimalSolution.stateTrajectory_;
            outputPrimalSolution.inputTrajectory_ = outputPrimalSolution.inputTrajectory_;

            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * Forward integrate the system dynamics with the controller in inputPrimalSolution. In general, it uses the given
     * control policies and the initial state, to integrate the system dynamics in the time period [initTime, finalTime].
     * However, if inputPrimalSolution's controller does not cover the period [initTime, finalTime], it will use the
     * controller till the final time of the controller
     *
     * @param [in] inputPrimalSolution: Its controller will be used for rollout.
     * @param [out] outputPrimalSolution: The resulting PrimalSolution.
     * @return True if the rollout was successful.
     */
    bool rolloutInitialController(PrimalSolution_t &inputPrimalSolution, PrimalSolution_t &outputPrimalSolution)
    {
        if (inputPrimalSolution.controllerPtr_ == nullptr)
        {
            return false;
        }

        // const Scalar finalTime = std::max(initTime_, std::min(inputPrimalSolution.controllerPtr_->timeStamp_.back(), finalTime_));

        if (initTime_ < finalTime_)
        {
            outputPrimalSolution.controllerPtr_ = inputPrimalSolution.controllerPtr_;
            rolloutTrajectory(*rolloutPtr_, initTime_, initState_, finalTime_, outputPrimalSolution);
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * Initializes the nominal dual solutions based on the optimized ones and nominal primal solution.
     * Moreover, it updates ProblemMetrics.
     */
    void initializeDualSolutionAndMetrics()
    {
        // initialize dual solution
        initializeDualSolution(optimalControlProblem_, nominalPrimalData_.primalSolution, optimizedDualSolution_,
                               nominalDualData_.dualSolution);

        computeRolloutMetrics(optimalControlProblem_, nominalPrimalData_.primalSolution, nominalDualData_.dualSolution,
                              nominalPrimalData_.problemMetrics);

        // update dual
        //  totalDualSolutionTimer_.startTimer();
        //  ocs2::updateDualSolution(optimalControlProblemStock_[0], nominalPrimalData_.primalSolution, nominalPrimalData_.problemMetrics,
        //  nominalDualData_.dualSolution);
        //  totalDualSolutionTimer_.endTimer();

        // calculates rollout merit
        performanceIndex_ = computeRolloutPerformanceIndex(nominalPrimalData_.primalSolution.timeTrajectory_, nominalPrimalData_.problemMetrics);
        performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
    }

    /**
     * Calculates the merit function based on the performance index .
     *
     * @param [in] performanceIndex: The performance index which includes the uninitialized merit, cost, and ISEs of constraints.
     * @return The merit function
     */
    Scalar calculateRolloutMerit(const PerformanceIndex_t &performanceIndex) const
    {
        // cost
        Scalar merit = performanceIndex.cost;
        // state/state-input equality Lagrangian
        merit += performanceIndex.equalityLagrangian;
        // state/state-input inequality Lagrangian
        merit += performanceIndex.inequalityLagrangian;

        return merit;
    }

    /**
     * Approximates the nonlinear problem as a linear-quadratic problem around the
     * nominal state and control trajectories. This method updates the following
     * variables:
     * 	- linearized system model and constraints
     * 	- quadratized cost function
     * 	- as well as the constrained coefficients of
     * 		- linearized system model
     * 		- quadratized intermediate cost function
     * 		- quadratized final cost
     */
    void approximateOptimalControlProblem()
    {
        approximateIntermediateLQ(nominalDualData_.dualSolution, nominalPrimalData_);

        /*
         * compute the Heuristics function at the final time. Also call shiftHessian on the Heuristics 2nd order derivative.
         */
        ModelData_t &modelData = nominalPrimalData_.modelDataFinalTime;
        const Scalar &time = nominalPrimalData_.primalSolution.timeTrajectory_.back();
        const StateVector_t &state = nominalPrimalData_.primalSolution.stateTrajectory_.back();
        const FinalMultiplierCollection_t &multiplier = nominalDualData_.dualSolution.final;
        modelData = approximator_.approximateFinalLQ(optimalControlProblem_, time, state, multiplier);

        // shift Hessian for final time
        if (ddpSettings_.strategy_ == SearchStrategyType::LINE_SEARCH)
        {
            shiftHessian(ddpSettings_.lineSearch_.hessianCorrectionStrategy, modelData.cost.dfdxx,
                         ddpSettings_.lineSearch_.hessianCorrectionMultiple);
        }
    }

    /**
     * Calculates an LQ approximate of the optimal control problem for the nodes.
     *
     * @param [in] dualSolution: The dual solution
     * @param [in,out] primalData: The primal Data
     */
    void approximateIntermediateLQ(const DualSolution_t &dualSolution, PrimalDataContainer_t &primalData)
    {
        // create alias
        const TimeTrajectory_t &timeTrajectory = primalData.primalSolution.timeTrajectory_;
        const StateTrajectory_t &stateTrajectory = primalData.primalSolution.stateTrajectory_;
        const InputTrajectory_t &inputTrajectory = primalData.primalSolution.inputTrajectory_;
        const IntermediateMultiplierTrajectory_t &multiplierTrajectory = dualSolution.intermediates;
        ModelDataTrajectory_t &modelDataTrajectory = primalData.modelDataTrajectory;

        for (size_t timeIndex = 0; timeIndex < PredictLength; ++timeIndex)
        {
            // approximate continuous LQ for the given time index
            ModelData_t continuousTimeModelData = approximator_.approximateIntermediateLQ(optimalControlProblem_, timeTrajectory[timeIndex], stateTrajectory[timeIndex],
                                                                                          inputTrajectory[timeIndex], multiplierTrajectory[timeIndex]);

            // TO DO:checking the numerical properties

            // discretize LQ problem
            const Scalar timeStep = (timeIndex + 1 < timeTrajectory.size()) ? (timeTrajectory[timeIndex + 1] - timeTrajectory[timeIndex]) : 0.0;
            if (!numerics::almost_eq(timeStep, 0.0))
            {
                discreteLQWorker(*optimalControlProblem_.dynamicsPtr, timeTrajectory[timeIndex], stateTrajectory[timeIndex],
                                 inputTrajectory[timeIndex], timeStep, continuousTimeModelData, modelDataTrajectory[timeIndex]);
            }
            else
            {
                modelDataTrajectory[timeIndex] = continuousTimeModelData;
            }
        };
    }

    /**
     * Calculates the discrete-time LQ approximation from the continuous-time LQ approximation.
     *
     * @param [in] system: system dynamic.
     * @param [in] time: time t_k.
     * @param [in] state: state x_k.
     * @param [in] input: input u_k.
     * @param [in] timeStep: Time step between the x_{k} and x_{k+1}.
     * @param [in] continuousTimeModelData: continuous time model data.
     * @param [out] modelData: Discretized mode data.
     */
    void discreteLQWorker(SystemDynamicsBase<Scalar, XDimisions, UDimisions> &system, Scalar time, const StateVector_t &state, const InputVector_t &input, Scalar timeStep,
                          const ModelData_t &continuousTimeModelData, ModelData_t &modelData)
    {
        modelData.time = continuousTimeModelData.time;

        // linearize system dynamics
        modelData.dynamics = discretizer_.sensitivityDiscretize(system, time, state, input, timeStep);
        modelData.dynamics.f.setZero(); // why?

        // quadratic approximation to the cost function
        modelData.cost = continuousTimeModelData.cost * timeStep;
    }

    /**
     * Computes the Hessian of Hamiltonian based on the search strategy and algorithm.
     *
     * @param [in] modelData: The model data.
     * @param [in] Sm: The Riccati matrix.
     * @return The Hessian matrix of the Hamiltonian.
     */
    HmMatrix_t computeHamiltonianHessian(const ModelData_t &modelData, const SmMatrix_t &Sm) const
    {
        const Matrix<Scalar, UDimisions, XDimisions> BmTransSm = modelData.dynamics.dfdu.transpose() * Sm;
        HmMatrix_t Hm = modelData.cost.dfduu;
        Hm += BmTransSm * modelData.dynamics.dfdu;
        return searchStrategyPtr_->augmentHamiltonianHessian(modelData, Hm);
    }

    /**
     *
     * @param [in] Hm: inv(Hm) defines the oblique projection for state-input equality constraints.
     * @param [out] constraintRangeProjector: The projection matrix to the constrained subspace.
     * @param [out] constraintNullProjector: The projection matrix to the null space of constrained.
     */
    void computeProjections(const HmMatrix_t &Hm, HmMatrix_t &constraintNullProjector) const
    {
        // UUT decomposition of inv(Hm)
        HmMatrix_t HmInvUmUmT;

        LinearAlgebra::computeInverseMatrixUUT(Hm, HmInvUmUmT);

        // compute DmDagger, DmDaggerTHmDmDaggerUUT, HmInverseConstrainedLowRank
        constraintNullProjector = HmInvUmUmT;
    }

    /**
     * Projects the unconstrained LQ coefficients to constrained ones.
     *
     * @param [in] modelData: The model data.
     * @param [in] constraintRangeProjector: The projection matrix to the constrained subspace.
     * @param [in] constraintNullProjector: The projection matrix to the null space of constrained.
     * @param [out] projectedModelData: The projected model data.
     */
    void projectLQ(const ModelData_t &modelData, const HmMatrix_t &constraintNullProjector, ModelData_t &projectedModelData) const
    {
        // dimensions and time
        projectedModelData.time = modelData.time;

        // Change of variables u = Pu * tilde{u}
        // Pu = constraintNullProjector;

        // dynamics
        projectedModelData.dynamics = modelData.dynamics;
        changeOfInputVariables(projectedModelData.dynamics, constraintNullProjector);

        // cost
        projectedModelData.cost = modelData.cost;
        changeOfInputVariables(projectedModelData.cost, constraintNullProjector);
    }

    /**
     * Takes the following steps: (1) Computes the Hessian of the Hamiltonian (i.e., Hm) (2) Based on Hm, it calculates
     * the range space and the null space projections of the input-state equality constraints. (3) Based on these two
     * projections, defines the projected LQ model. (4) Finally, defines the Riccati equation modifiers based on the
     * search strategy.
     *
     * @param [in] modelData: The model data.
     * @param [in] Sm: The Riccati matrix.
     * @param [out] projectedModelData: The projected model data.
     * @param [out] riccatiModification: The Riccati equation modifier.
     */
    void computeProjectionAndRiccatiModification(const ModelData_t &modelData, const SmMatrix_t &Sm, ModelData_t &projectedModelData,
                                                 RiccatiModification_t &riccatiModification) const
    {
        // compute the Hamiltonian's Hessian
        riccatiModification.time_ = modelData.time;
        riccatiModification.hamiltonianHessian_ = computeHamiltonianHessian(modelData, Sm);

        // compute projectors
        computeProjections(riccatiModification.hamiltonianHessian_, riccatiModification.constraintNullProjector_);

        // project LQ
        projectLQ(modelData, riccatiModification.constraintNullProjector_, projectedModelData);

        // compute deltaQm, deltaGv, deltaGm
        searchStrategyPtr_->computeRiccatiModification(projectedModelData, riccatiModification.deltaQm_);
    }

    /**
     * Solves Riccati equations.

        // dynamics bias
        projectedModelData.dynamicsBias = modelDat
     *
     * @param [in] finalValueFunction: The final value of Sm (dfdxx), Sv (dfdx), s (f), for Riccati equation.
     * @return average time step
     */
    Scalar solveSequentialRiccatiEquations(const ValueFunctionQuadraticApproximation_t &finalValueFunction)
    {
        const ModelData_t &finalModelData = nominalPrimalData_.modelDataTrajectory.back();
        RiccatiModification_t &finalRiccatiModification = nominalDualData_.riccatiModificationTrajectory.back();
        ModelData_t &finalProjectedModelData = nominalDualData_.projectedModelDataTrajectory.back();
        LvVector_t &finalProjectedLvFinal = projectedLvTrajectoryStock_.back();
        KmMatrix_t &finalProjectedKmFinal = projectedKmTrajectoryStock_.back();

        Matrix<Scalar, XDimisions, XDimisions> SmDummy;
        SmDummy.setZero();

        computeProjectionAndRiccatiModification(finalModelData, SmDummy, finalProjectedModelData, finalRiccatiModification);

        // projected feedforward
        finalProjectedLvFinal = -finalProjectedModelData.cost.dfdu;
        finalProjectedLvFinal -= finalProjectedModelData.dynamics.dfdu.transpose() * finalValueFunction.dfdx;

        // projected feedback
        finalProjectedKmFinal = -finalProjectedModelData.cost.dfdux;
        finalProjectedKmFinal -= finalProjectedModelData.dynamics.dfdu.transpose() * finalValueFunction.dfdxx;

        return solveSequentialRiccatiEquationsImpl(finalValueFunction);
    }

    /**
     * The implementation for solving Riccati equations for all the partitions.
     *
     * @param [in] finalValueFunction The final Sm(dfdxx), Sv(dfdx), s(f), for Riccati equation.
     * @return average time step
     */
    Scalar solveSequentialRiccatiEquationsImpl(const ValueFunctionQuadraticApproximation_t &finalValueFunction)
    {
        // the last index of the partition is excluded, namely [first, last), so the value function approximation of the end point of the end
        // partition is filled manually.
        // For other partitions except the last one, the end points are filled in the solving stage of the next partition. For example,
        // [first1,last1), [first2(last1), last2).
        nominalDualData_.valueFunctionTrajectory.back() = finalValueFunction;

        // solve it sequentially for the first iteration
        const std::pair<int, int> partitionInterval{0, PredictLength - 1};
        riccatiEquationsWorker(partitionInterval, finalValueFunction);

        // average time step
        return (finalTime_ - initTime_) / PredictLength;
    }

    /**
     * Solves a Riccati equations and type_1 constraints error correction compensation for the partition in the given index.
     *
     * @param [in] workerIndex: Current worker index
     * @param [in] partitionInterval: Current active interval
     * @param [in] finalValueFunction The final Sm(dfdxx), Sv(dfdx), s(f), for Riccati equation.
     */
    void riccatiEquationsWorker(const std::pair<int, int> &partitionInterval, const ValueFunctionQuadraticApproximation_t &finalValueFunction)
    {
        // final temporal values. Used to store pre-jump value
        ValueFunctionQuadraticApproximation_t finalValueTemp = finalValueFunction;

        /*
         * solving the Riccati equations
         */
        const ValueFunctionQuadraticApproximation_t *valueFunctionNext = &finalValueTemp;

        int curIndex = partitionInterval.second - 1;
        const int stopIndex = partitionInterval.first;
        while (curIndex >= stopIndex)
        {
            LvVector_t &curProjectedLv = projectedLvTrajectoryStock_[curIndex];
            KmMatrix_t &curProjectedKm = projectedKmTrajectoryStock_[curIndex];
            ModelData_t &curProjectedModelData = nominalDualData_.projectedModelDataTrajectory[curIndex];
            RiccatiModification_t &curRiccatiModification = nominalDualData_.riccatiModificationTrajectory[curIndex];
            const ModelData_t &curModelData = nominalPrimalData_.modelDataTrajectory[curIndex];

            SmMatrix_t &curSm = nominalDualData_.valueFunctionTrajectory[curIndex].dfdxx;
            SvVector_t &curSv = nominalDualData_.valueFunctionTrajectory[curIndex].dfdx;
            Scalar &curs = nominalDualData_.valueFunctionTrajectory[curIndex].f;

            computeProjectionAndRiccatiModification(curModelData, valueFunctionNext->dfdxx, curProjectedModelData, curRiccatiModification);

            riccatiEquationsSolver_.computeMap(curProjectedModelData, curRiccatiModification, valueFunctionNext->dfdxx,
                                               valueFunctionNext->dfdx, valueFunctionNext->f, curProjectedKm, curProjectedLv, curSm,
                                               curSv, curs);
            valueFunctionNext = &(nominalDualData_.valueFunctionTrajectory[curIndex]);

            --curIndex;
        } // while
    }

    /**
     * Calculates the controller. This method uses the following variables. The method modifies unoptimizedController_.
     */
    void calculateController()
    {
        const size_t N = nominalPrimalData_.primalSolution.timeTrajectory_.size();

        unoptimizedController_.timeStamp_ = nominalPrimalData_.primalSolution.timeTrajectory_;

        size_t timeIndex = 0;
        // get next time index (atomic)
        while (timeIndex < N)
        {
            calculateControllerWorker(timeIndex, nominalPrimalData_, nominalDualData_, unoptimizedController_);
            timeIndex++;
        }

        // Since the controller for the last timestamp is invalid, if the last time is not the event time, use the control policy of the second to
        // last time for the last time
        // finalTimeIsNotAnEvent && there are at least two time stamps
        if (unoptimizedController_.size() >= 2)
        {
            const size_t secondToLastIndex = unoptimizedController_.size() - 2u;
            unoptimizedController_.gainArray_.back() = unoptimizedController_.gainArray_[secondToLastIndex];
            unoptimizedController_.biasArray_.back() = unoptimizedController_.biasArray_[secondToLastIndex];
            unoptimizedController_.deltaBiasArray_.back() = unoptimizedController_.deltaBiasArray_[secondToLastIndex];
        }
    }

    /**
     * Calculate controller for the timeIndex by using primal and dual and write the result back to dstController
     *
     * @param [in] timeIndex: The current time index
     * @param [in] primalData: Primal data used to calculate controller
     * @param [in] dualData: Dual data used to calculate controller
     * @param [out] dstController: The destination controller
     */
    void calculateControllerWorker(size_t timeIndex, const PrimalDataContainer_t &primalData, const DualDataContainer_t &dualData,
                                   LinearController_t &dstController)
    {
        const StateVector_t &nominalState = primalData.primalSolution.stateTrajectory_[timeIndex];
        const InputVector_t &nominalInput = primalData.primalSolution.inputTrajectory_[timeIndex];

        const Matrix<Scalar, UDimisions, UDimisions> &Qu = dualData.riccatiModificationTrajectory[timeIndex].constraintNullProjector_;

        // feedback gains
        dstController.gainArray_[timeIndex] = Qu * projectedKmTrajectoryStock_[timeIndex];

        // bias input
        dstController.biasArray_[timeIndex] = nominalInput;
        dstController.biasArray_[timeIndex] -= dstController.gainArray_[timeIndex] * nominalState;
        dstController.deltaBiasArray_[timeIndex] = Qu * projectedLvTrajectoryStock_[timeIndex];
    }

    /** Based on the current LQ solution updates the optimized primal and dual solutions. */
    void takePrimalDualStep(Scalar lqModelExpectedCost)
    {
        // update primal: run search strategy and find the optimal stepLength
        Scalar avgTimeStep = 0;
        SearchStrategySolutionRef_t solution(avgTimeStep, optimizedDualSolution_, optimizedPrimalSolution_, optimizedProblemMetrics_,
                                             performanceIndex_);
        const bool success = searchStrategyPtr_->run({initTime_, finalTime_}, initState_, lqModelExpectedCost, unoptimizedController_,
                                                     nominalDualData_.dualSolution, solution);

        if (success)
        {
            avgTimeStepFP_ = 0.9 * avgTimeStepFP_ + 0.1 * avgTimeStep;
        }

        // update dual
        if (success)
        {
            DualSolutionRef_t DualSolutionRef = optimizedDualSolution_;
            updateDualSolution(optimalControlProblem_, optimizedPrimalSolution_, optimizedProblemMetrics_, DualSolutionRef);
            performanceIndex_ = computeRolloutPerformanceIndex(optimizedPrimalSolution_.timeTrajectory_, optimizedProblemMetrics_);
            performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
        }

        // if failed, use nominal and to keep the consistency of cached data, all cache should be left untouched
        if (!success)
        {
            optimizedDualSolution_ = nominalDualData_.dualSolution;
            optimizedPrimalSolution_ = nominalPrimalData_.primalSolution;
            optimizedProblemMetrics_ = nominalPrimalData_.problemMetrics;
            performanceIndex_ = performanceIndexLast_;
        }
    }
    /**
     * Computes cost, soft constraints and constraints values of each point in the the primalSolution rollout.
     *
     * @param [in] problem: A reference to the optimal control problem.
     * @param [in] primalSolution: The primal solution.
     * @param [in] dualSolution: Const reference view to the dual solution
     * @param [out] problemMetrics: The cost, soft constraints and constraints values of the rollout.
     */
    void computeRolloutMetrics(OptimalControlProblem_t &problem, const PrimalSolution_t &primalSolution,
                               DualSolution_t &dualSolution, ProblemMetrics_t &problemMetrics)
    {
        const std::array<Scalar, PredictLength> &tTrajectory = primalSolution.timeTrajectory_;
        const std::array<Vector<Scalar, XDimisions>, PredictLength> &xTrajectory = primalSolution.stateTrajectory_;
        const std::array<Vector<Scalar, UDimisions>, PredictLength> &uTrajectory = primalSolution.inputTrajectory_;

        for (size_t k = 0; k < PredictLength; k++)
        {
            // intermediate time cost and constraints
            problemMetrics.intermediates[k] =
                approximator_.computeIntermediateMetrics(problem, tTrajectory[k], xTrajectory[k], uTrajectory[k], dualSolution.intermediates[k]);
        }

        // final time cost and constraints
        problemMetrics.final = approximator_.computeFinalMetrics(problem, tTrajectory.back(), xTrajectory.back(), dualSolution.final);
    }

private:
    OptimalControlProblem_t optimalControlProblem_;

    DDPSettings<Scalar> ddpSettings_;

    // time and state
    Scalar initTime_ = 0.0;
    Scalar finalTime_ = 0.0;
    StateVector_t initState_;

    // nominal data
    PrimalDataContainer_t nominalPrimalData_;
    DualDataContainer_t nominalDualData_;

    // controller that is calculated directly from dual solution. It is unoptimized because it haven't gone through searching.
    LinearController_t unoptimizedController_;

    SearchStrategyBase_t *searchStrategyPtr_;

    // reference trajectory
    TimeTrajectory_t timeTrajectory_;
    StateTrajectory_t stateTrajectory_;
    InputTrajectory_t inputTrajectory_;

    // optimized data
    PrimalSolution_t optimizedPrimalSolution_;
    DualSolution_t optimizedDualSolution_;
    ProblemMetrics_t optimizedProblemMetrics_;

    // performance index
    PerformanceIndex_t performanceIndex_;
    PerformanceIndex_t performanceIndexLast_;

    // roll out
    RolloutBase_t *rolloutPtr_;
    TimeTriggeredRollout_t rollout_;

    // linear approximator
    LinearQuadraticApproximator_t approximator_;

    // Discretizer
    EK2DynamicsDiscretizer_t discretizer_;

    // Discrete time riccati equation solver
    DiscreteTimeRiccatiEquations_t riccatiEquationsSolver_;
    std::array<KmMatrix_t, PredictLength> projectedKmTrajectoryStock_; // projected feedback
    std::array<LvVector_t, PredictLength> projectedLvTrajectoryStock_; // projected feedforward

    // forward pass and backward pass average time step
    Scalar avgTimeStepFP_ = 0.0;
    Scalar avgTimeStepBP_ = 0.0;

    size_t totalNumIterations_{0};

    Scalar predictTimeStep_;
};