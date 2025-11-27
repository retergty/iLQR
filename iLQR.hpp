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
#include "TrapezoidalIntegration.hpp"
#include "Initializer.hpp"
#include "InitializerRollout.hpp"

template <typename Scalar, int XDimisions, int UDimisions, size_t PredictLength,
          int StateEqConstrains = 0, int StateIneqConstrains = 0, int StateInputEqConstrains = 0, int StateInputIneqConstrains = 0,
          int FinalStateEqConstrains = 0, int FinalStateIneqConstrains = 0>
class iLQR
{
public:
    using OptimalControlProblem_t = OptimalControlProblem<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using ControlledSystemBase_t = ControlledSystemBase<Scalar, XDimisions, UDimisions>;
    using SystemDynamicsBase_t = SystemDynamicsBase<Scalar, XDimisions, UDimisions>;

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
    using RolloutBase_t = RolloutBase<Scalar, XDimisions, UDimisions>;
    using InitializerRollout_t = InitializerRollout<Scalar, XDimisions, UDimisions>;
    using Initializer_t = Initializer<Scalar, XDimisions, UDimisions>;

    using TimeTriggeredRollout_t = TimeTriggeredRollout<Scalar, XDimisions, UDimisions>;
    using RolloutTrajectoryPointer_t = RolloutTrajectoryPointer<Scalar, XDimisions, UDimisions>;

    using SearchStrategySolution_t = SearchStrategySolution<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using SearchStrategySolutionRef_t = SearchStrategySolutionRef<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;

    using LinearController_t = LinearController<Scalar, XDimisions, UDimisions, PredictLength + 1>;
    using SearchStrategyBase_t = SearchStrategyBase<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using LineSearchStrategy_t = LineSearchStrategy<Scalar, XDimisions, UDimisions, PredictLength, StateEqConstrains, StateIneqConstrains, StateInputEqConstrains, StateInputIneqConstrains, FinalStateEqConstrains, FinalStateIneqConstrains>;
    using RiccatiModification_t = RiccatiModification<Scalar, XDimisions, UDimisions>;

    using TimeTrajectory_t = std::array<Scalar, PredictLength + 1>;
    using StateTrajectory_t = std::array<Vector<Scalar, XDimisions>, PredictLength + 1>;
    using InputTrajectory_t = std::array<Vector<Scalar, UDimisions>, PredictLength + 1>;
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
    using IntermediatePerformanceIndexTrajectory_t = std::array<PerformanceIndex_t, PredictLength>;

    using EK2DynamicsDiscretizer_t = EK2DynamicsDiscretizer<Scalar, XDimisions, UDimisions>;

    using ValueFunctionQuadraticApproximation_t = ScalarFunctionQuadraticApproximation<Scalar, XDimisions, UDimisions>;
    using DiscreteTimeRiccatiEquations_t = DiscreteTimeRiccatiEquations<Scalar, XDimisions, UDimisions>;

    iLQR(SystemDynamicsBase_t *systemPtr, Initializer_t *initializer) : rollout_(systemPtr, ddpSettings_.timeStep_),
                                                                        initializerRollout_(*initializer, ddpSettings_.timeStep_),
                                                                        lineSearchStrategy_(*this)
    {
        optimalControlProblem_.dynamicsPtr = systemPtr;
        // set zero solution
        optimizedPrimalSolution_.clear();
        optimizedDualSolution_.clear();
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
        lastFinalTime_ = finalTime_;
        finalTime_ = initTime + ddpSettings_.timeStep_ * (PredictLength);
        const size_t initIteration = totalNumIterations_;

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
            isConverged = lineSearchStrategy_.checkConvergence(
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
        int numSteps = 0;
        if (lastFinalTime_ > initTime_)
        {
            numSteps = rolloutInitialController(optimizedPrimalSolution_, nominalPrimalData_.primalSolution);
        }
        // past too long before last run.
        rolloutInitializer(nominalPrimalData_.primalSolution, numSteps);
        return true;
    }

    /**
     * Forward integrate the system dynamics with the controller in inputPrimalSolution. In general, it uses the given
     * control policies and the initial state, to integrate the system dynamics in the time period [initTime, finalTime].
     * However, if inputPrimalSolution's controller does not cover the period [initTime, finalTime], it will use the
     * controller till the final time of the controller
     *
     * @param [in] inputPrimalSolution: Its controller will be used for rollout.
     * @param [out] outputPrimalSolution: The resulting PrimalSolution.
     * @return number of steps are covered.
     */
    int rolloutInitialController(PrimalSolution_t &inputPrimalSolution, PrimalSolution_t &outputPrimalSolution)
    {
        // Ensure that finalTime is included by adding a fraction of dt such that: N * dt <= finalTime < (N + 1) * dt.
        Scalar finalTimeLocal = std::min(lastFinalTime_, finalTime_) + static_cast<Scalar>(0.01) * ddpSettings_.timeStep_;
        int numSteps = (finalTimeLocal - initTime_) / PredictLength;

        outputPrimalSolution.controller_ = inputPrimalSolution.controller_;
        rolloutTrajectory(rollout_, initTime_, initState_, lastFinalTime_, outputPrimalSolution);

        return numSteps;
    }

    /**
     * It will check the content of the primalSolution, and if its final time is smaller than the current solver finalTime_,
     * it will concatenate it with the result of Initializer.
     */
    void rolloutInitializer(PrimalSolution_t &primalSolution, int numSteps)
    {
        RolloutTrajectoryPointer_t rolloutTrajectortPtr(primalSolution.timeTrajectory_.data() + numSteps, primalSolution.stateTrajectory_.data() + numSteps + 1, primalSolution.inputTrajectory_.data() + numSteps, PredictLength - numSteps + 1);
        Scalar initTime = lastFinalTime_;
        const StateVector_t &initState = primalSolution.stateTrajectory_[numSteps];

        initializerRollout_.run(initTime, initState, finalTime_, nullptr, rolloutTrajectortPtr);
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

        // calculates rollout merit
        performanceIndex_ = computeRolloutPerformanceIndex(nominalPrimalData_.primalSolution.timeTrajectory_, nominalPrimalData_.problemMetrics);
        performanceIndex_.merit = calculateRolloutMerit(performanceIndex_);
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

            // // discretize LQ problem
            // const Scalar timeStep = (timeIndex + 1 < timeTrajectory.size()) ? (timeTrajectory[timeIndex + 1] - timeTrajectory[timeIndex]) : 0.0;
            // if (!numerics::almost_eq(timeStep, 0.0))
            // {
            //     discreteLQWorker(*optimalControlProblem_.dynamicsPtr, timeTrajectory[timeIndex], stateTrajectory[timeIndex],
            //                      inputTrajectory[timeIndex], timeStep, continuousTimeModelData, modelDataTrajectory[timeIndex]);
            // }
            // else
            // {
            //     modelDataTrajectory[timeIndex] = continuousTimeModelData;
            // }

            const Scalar timeStep = timeTrajectory[timeIndex + 1] - timeTrajectory[timeIndex];
            discreteLQWorker(*optimalControlProblem_.dynamicsPtr, timeTrajectory[timeIndex], stateTrajectory[timeIndex],
                             inputTrajectory[timeIndex], timeStep, continuousTimeModelData, modelDataTrajectory[timeIndex]);
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
        return lineSearchStrategy_.augmentHamiltonianHessian(modelData, Hm);
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
        lineSearchStrategy_.computeRiccatiModification(projectedModelData, riccatiModification.deltaQm_);
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
        const ModelData_t &finalModelData = nominalPrimalData_.modelDataFinalTime;
        RiccatiModification_t &finalRiccatiModification = nominalDualData_.riccatiModificationTrajectory.back();
        ModelData_t &finalProjectedModelData = nominalDualData_.projectedModelDataTrajectory.back();
        LvVector_t &finalProjectedLvFinal = projectedLvTrajectoryStock_.back();
        KmMatrix_t &finalProjectedKmFinal = projectedKmTrajectoryStock_.back();

        SmMatrix_t SmDummy;
        SmDummy.setZero();

        computeProjectionAndRiccatiModification(finalModelData, SmDummy, finalProjectedModelData, finalRiccatiModification);

        // projected feedforward
        finalProjectedLvFinal = -finalProjectedModelData.cost.dfdu;
        // last
        // finalProjectedLvFinal -= finalProjectedModelData.dynamics.dfdu.transpose() * finalValueFunction.dfdx;

        // projected feedback
        finalProjectedKmFinal = -finalProjectedModelData.cost.dfdux;
        // finalProjectedKmFinal -= finalProjectedModelData.dynamics.dfdu.transpose() * finalValueFunction.dfdxx;

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
        nominalDualData_.valueFunctionTrajectory.back() = finalValueFunction;

        riccatiEquationsWorker(finalValueFunction);

        // average time step
        return (finalTime_ - initTime_) / PredictLength;
    }

    /**
     * Solves a Riccati equations and type_1 constraints error correction compensation for the partition in the given index.
     *
     * @param [in] workerIndex: Current worker index
     * @param [in] finalValueFunction The final Sm(dfdxx), Sv(dfdx), s(f), for Riccati equation.
     */
    void riccatiEquationsWorker(const ValueFunctionQuadraticApproximation_t &finalValueFunction)
    {
        /*
         * solving the Riccati equations
         */
        const ValueFunctionQuadraticApproximation_t *valueFunctionNext = &finalValueFunction;

        int curIndex = PredictLength - 1;
        const int stopIndex = 0;
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
        unoptimizedController_.timeStamp_ = nominalPrimalData_.primalSolution.timeTrajectory_;

        for (size_t timeIndex = 0; timeIndex < PredictLength; ++timeIndex)
        {
            calculateControllerWorker(timeIndex, nominalPrimalData_, nominalDualData_, unoptimizedController_);
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
        const bool success = lineSearchStrategy_.run({initTime_, finalTime_}, initState_, lqModelExpectedCost, unoptimizedController_,
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

public:
    /**
     * Calculates the merit function based on the performance index .
     *
     * @param [in] performanceIndex: The performance index which includes the uninitialized merit, cost, and ISEs of constraints.
     * @return The merit function
     */
    static Scalar calculateRolloutMerit(const PerformanceIndex_t &performanceIndex)
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
     * Computes cost, soft constraints and constraints values of each point in the the primalSolution rollout.
     *
     * @param [in] problem: A reference to the optimal control problem.
     * @param [in] primalSolution: The primal solution.
     * @param [in] dualSolution: Const reference view to the dual solution
     * @param [out] problemMetrics: The cost, soft constraints and constraints values of the rollout.
     */
    static void computeRolloutMetrics(OptimalControlProblem_t &problem, const PrimalSolution_t &primalSolution,
                                      DualSolution_t &dualSolution, ProblemMetrics_t &problemMetrics)
    {
        const TimeTrajectory_t &tTrajectory = primalSolution.timeTrajectory_;
        const StateTrajectory_t &xTrajectory = primalSolution.stateTrajectory_;
        const InputTrajectory_t &uTrajectory = primalSolution.inputTrajectory_;

        for (size_t k = 0; k < PredictLength; k++)
        {
            // intermediate time cost and constraints
            problemMetrics.intermediates[k] =
                approximator_.computeIntermediateMetrics(problem, tTrajectory[k], xTrajectory[k], uTrajectory[k], dualSolution.intermediates[k]);
        }

        // final time cost and constraints
        problemMetrics.final = approximator_.computeFinalMetrics(problem, tTrajectory.back(), xTrajectory.back(), dualSolution.final);
    }

    /**
     * Forward integrate the system dynamics with given controller. It uses the given control policies and initial state,
     * to integrate the system dynamics in time period [initTime, finalTime].
     *
     * @param [in] rollout: A reference to the rollout class.
     * @param [in] initTime: The initial time.
     * @param [in] initState: The initial state.
     * @param [in] finalTime: The final time.
     * @param [in, out] primalSolution: The resulting primal solution. Make sure that primalSolution::controllerPtr is set since
     *                                  the rollout is performed based on the controller stored in primalSolution. Moreover,
     *                                  except for StateTriggeredRollout, one should also set primalSolution::modeSchedule.
     *
     * @return average time step.
     */
    static Scalar rolloutTrajectory(RolloutBase_t &rollout,
                                    Scalar initTime, const StateVector_t &initState, Scalar finalTime,
                                    PrimalSolution_t &primalSolution)
    {
        RolloutTrajectoryPointer_t rolloutTrajectortPtr(primalSolution.timeTrajectory_.data(), primalSolution.stateTrajectory_.data(), primalSolution.inputTrajectory_.data(), PredictLength + 1);
        rollout.run(initTime, initState, finalTime, &primalSolution.controller_, rolloutTrajectortPtr);
        // average time step
        return (finalTime - initTime) / static_cast<Scalar>(PredictLength);
    }

    /**
     * Calculates the PerformanceIndex associated to the given ProblemMetrics.
     *
     * @param [in] timeTrajectory: Time stamp of the rollout.
     * @param [in] problemMetrics: The cost, soft constraints and constraints values of the rollout.
     *
     * @return The PerformanceIndex of the trajectory.
     */
    static PerformanceIndex_t computeRolloutPerformanceIndex(
        const TimeTrajectory_t &timeTrajectory,
        const ProblemMetrics_t &problemMetrics)
    {
        // Final
        PerformanceIndex_t finalperformanceIndex = toPerformanceIndex(problemMetrics.final);
        // Intermediates
        IntermediatePerformanceIndexTrajectory_t performanceIndexTrajectory;

        // std::for_each(problemMetrics.intermediates.cbegin(), problemMetrics.intermediates.cend(),
        //   [&performanceIndexTrajectory](const Metrics& m) { performanceIndexTrajectory.push_back(toPerformanceIndex(m)); });

        for (size_t i = 0; i < performanceIndexTrajectory.size(); ++i)
        {
            performanceIndexTrajectory[i] = toPerformanceIndex(problemMetrics.intermediates[i]);
        }

        // Intermediates
        return trapezoidalIntegration(timeTrajectory, performanceIndexTrajectory, finalperformanceIndex);
    }

    /**
     * Outputs a controller with the same time stamp and gains as unoptimizedController. However, bias is incremented based on:
     * biasArray = unoptimizedController.biasArray + stepLength * unoptimizedController.deltaBiasArray
     */
    static void incrementController(Scalar stepLength, const LinearController_t &unoptimizedController, LinearController_t &controller)
    {
        controller.timeStamp_ = unoptimizedController.timeStamp_;
        controller.gainArray_ = unoptimizedController.gainArray_;
        for (size_t k = 0; k < unoptimizedController.size(); k++)
        {
            controller.biasArray_[k] = unoptimizedController.biasArray_[k] + stepLength * unoptimizedController.deltaBiasArray_[k];
        }
    }

    /**
     * biasArray = unoptimizedController.biasArray + stepLength * unoptimizedController.deltaBiasArray
     */
    static void changeControllerStepLength(Scalar stepLength, const LinearController_t &unoptimizedController, LinearController_t &controller)
    {
        for (size_t k = 0; k < unoptimizedController.size(); k++)
        {
            controller.biasArray_[k] = unoptimizedController.biasArray_[k] + stepLength * unoptimizedController.deltaBiasArray_[k];
        }
    }

    /**
     * Computes the integral of the squared (IS) norm of the controller update.
     *
     * @param [in] controller: Input controller.
     * @return The integral of the squared (IS) norm of the controller update.
     */
    static Scalar computeControllerUpdateIS(const LinearController_t &controller)
    {
        std::array<Scalar, controller.size()> biasArraySquaredNorm;

        for (size_t i = 0; i < controller.size(); ++i)
        {
            biasArraySquaredNorm[i] = controller.deltaBiasArray_[i].squaredNorm();
        }
        // integrates using the trapezoidal approximation method
        return trapezoidalIntegration(controller.timeStamp_, biasArraySquaredNorm, 0.0);
    }

    OptimalControlProblem_t optimalControlProblem_;

    // roll out
    TimeTriggeredRollout_t rollout_;

    // initializer
    InitializerRollout_t initializerRollout_;

private:
    // linear approximator
    static LinearQuadraticApproximator_t approximator_;

    constexpr static DDPSettings<Scalar> ddpSettings_{};

    // time and state
    Scalar initTime_{0.0};
    Scalar finalTime_{0.0};
    Scalar lastFinalTime_{0.0};

    StateVector_t initState_;

    // nominal data
    PrimalDataContainer_t nominalPrimalData_;
    DualDataContainer_t nominalDualData_;

    // controller that is calculated directly from dual solution. It is unoptimized because it haven't gone through searching.
    LinearController_t unoptimizedController_;

    LineSearchStrategy_t lineSearchStrategy_;

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

    // Discretizer
    EK2DynamicsDiscretizer_t discretizer_;

    // Discrete time riccati equation solver
    DiscreteTimeRiccatiEquations_t riccatiEquationsSolver_;
    std::array<KmMatrix_t, PredictLength + 1> projectedKmTrajectoryStock_; // projected feedback
    std::array<LvVector_t, PredictLength + 1> projectedLvTrajectoryStock_; // projected feedforward

    // forward pass and backward pass average time step
    Scalar avgTimeStepFP_ = 0.0;
    Scalar avgTimeStepBP_ = 0.0;

    size_t totalNumIterations_{0};
};