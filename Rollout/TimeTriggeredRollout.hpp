#pragma once
#include "RolloutBase.hpp"
#include "RungeKuttaDormandPrince5.hpp"

/**
 * This class is an interface class for forward rollout of the system dynamics.
 */
template <typename Scalar, int XDimisions, int UDimisions>
class TimeTriggeredRollout : public RolloutBase<Scalar, XDimisions, UDimisions>
{
public:
    using RolloutTrajectoryPointer_t = typename RolloutBase<Scalar, XDimisions, UDimisions>::RolloutTrajectoryPointer_t;
    /**
     * Constructor.
     *
     * @param [in] systemDynamics: The system dynamics for forward rollout.
     * @param [in] rolloutSettings: The rollout settings.
     */
    explicit TimeTriggeredRollout(ControlledSystemBase<Scalar, XDimisions, UDimisions>* systemDynamics, const Scalar timeStep)
        : systemDynamicsPtr_(systemDynamics)
    {
        this->rolloutSettings_.timeStep = timeStep;
    };

    ~TimeTriggeredRollout() override = default;

    /** Returns the underlying dynamics. */
    ControlledSystemBase<Scalar, XDimisions, UDimisions>* systemDynamicsPtr() { return systemDynamicsPtr_; }

    int run(const Scalar initTime, const Vector<Scalar, XDimisions>& initState, const Scalar finalTime, ControllerBase<Scalar, XDimisions, UDimisions>* controller,
        RolloutTrajectoryPointer_t& trajectory) override
    {
        assert(finalTime > initTime);

        // set controller
        systemDynamicsPtr_->setController(controller);

        Observer<Scalar, XDimisions> observer(trajectory.maxLength, trajectory.stateTrajectory, trajectory.timeTrajectory); // concatenate trajectory
        // integrate controlled system
        RK45Intergraor_.integrateConst(*systemDynamicsPtr_, observer, initState, initTime, finalTime, this->settings().timeStep);

        int RolloutIntegrateCount = observer.getCount();

        // compute control input trajectory and concatenate to inputTrajectory
        if (this->settings().reconstructInputTrajectory)
        {
            for (int i = 0; i < RolloutIntegrateCount; i++)
            {
                trajectory.inputTrajectory[i] = systemDynamicsPtr_->controllerPtr()->computeInput(trajectory.timeTrajectory[i], trajectory.stateTrajectory[i]);
            } // end of k_u loop
        }

        return RolloutIntegrateCount;
    }

    // Scalar run(const Scalar initTime, const Vector<Scalar, XDimisions>& initState, const int steps, ControllerBase<Scalar, XDimisions, UDimisions>& controller,
    //     std::array<Scalar, ArrayLen>& timeTrajectory, std::array<Vector<Scalar, XDimisions>, ArrayLen>& stateTrajectory, std::array<Vector<Scalar, UDimisions>, ArrayLen>& inputTrajectory) override
    // {
    //     Scalar finalTime = initTime + steps * this->settings().timeStep;
    //     // set controller
    //     systemDynamicsPtr_->setController(&controller);

    //     Observer<Scalar, XDimisions> observer(ArrayLen, stateTrajectory.data(), timeTrajectory.data()); // concatenate trajectory
    //     // integrate controlled system
    //     RK45Intergraor_.integrateConst(*systemDynamicsPtr_, observer, initState, initTime, finalTime, this->settings().timeStep);
    // }

private:
    ControlledSystemBase<Scalar, XDimisions, UDimisions>* systemDynamicsPtr_{ nullptr };

    RungeKuttaDormandPrince5<Scalar, XDimisions> RK45Intergraor_;
};
