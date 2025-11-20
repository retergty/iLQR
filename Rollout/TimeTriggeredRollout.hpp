#pragma once
#include "RolloutBase.hpp"

/**
 * This class is an interface class for forward rollout of the system dynamics.
 */
template <typename Scalar, int XDimisions, int UDimisions, size_t ArrayLen>
class TimeTriggeredRollout : public RolloutBase<Scalar, XDimisions, UDimisions, ArrayLen>
{
public:
    /**
     * Constructor.
     *
     * @param [in] systemDynamics: The system dynamics for forward rollout.
     * @param [in] rolloutSettings: The rollout settings.
     */
    explicit TimeTriggeredRollout(ControlledSystemBase<Scalar, XDimisions, UDimisions>* systemDynamics, IntegratorBase<Scalar, XDimisions>* dynamic_integrator, const RolloutSettings<Scalar>& rolloutSettings = RolloutSettings())
        : system_dynamics_ptr_(systemDynamics), dynamics_integrator_ptr_(dynamic_integrator), RolloutBase(rolloutSettings) {
    };

    ~TimeTriggeredRollout() override = default;
    TimeTriggeredRollout(const TimeTriggeredRollout&) = delete;
    TimeTriggeredRollout& operator=(const TimeTriggeredRollout&) = delete;

    /** Returns the underlying dynamics. */
    ControlledSystemBase* systemDynamicsPtr() { return systemDynamicsPtr_.get(); }

    int run(const Scalar initTime, const Vector<Scalar, XDimisions>& initState, const Scalar finalTime, ControllerBase<Scalar, XDimisions, UDimisions>* controller,
        std::array<Scalar, ArrayLen>& timeTrajectory, std::array<Vector<Scalar, XDimisions>, ArrayLen>& stateTrajectory, std::array<Vector<Scalar, UDimisions>, ArrayLen>& inputTrajectory) override
    {
        assert(finalTime > initTime);
        assert(controller != nullptr);

        // set controller
        system_dynamics_ptr_->setController(controller);

        Observer<Scalar, XDimisions> observer(ArrayLen, stateTrajectory.data(), timeTrajectory.data()); // concatenate trajectory
        // integrate controlled system
        dynamics_integrator_ptr_->integrateConst(*systemDynamicsPtr_, observer, initState, initTime, finalTime, this->settings().timeStep);

        int RolloutIntegrateCount = observer.getCount();

        // compute control input trajectory and concatenate to inputTrajectory
        if (this->settings().reconstructInputTrajectory)
        {
            for (int i = 0; i < RolloutIntegrateCount; i++)
            {
                inputTrajectory[i] = systemDynamicsPtr_->controllerPtr()->computeInput(timeTrajectory[i], stateTrajectory[i]);
            } // end of k_u loop
        }

        return RolloutIntegrateCount;
    }

private:
    ControlledSystemBase<Scalar, XDimisions, UDimisions>* system_dynamics_ptr_{ nullptr };

    IntegratorBase<Scalar, XDimisions>* dynamics_integrator_ptr_{ nullptr };
};
