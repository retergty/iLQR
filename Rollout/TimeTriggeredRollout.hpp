#pragma once
#include "RolloutBase.hpp"
#include "RungeKuttaDormandPrince5.hpp"

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
    explicit TimeTriggeredRollout(ControlledSystemBase<Scalar, XDimisions, UDimisions> *systemDynamics, const Scalar timeStep)
        : systemDynamicsPtr_(systemDynamics)
    {
        this->rolloutSettings_.timeStep = timeStep;
    };

    ~TimeTriggeredRollout() override = default;

    /** Returns the underlying dynamics. */
    ControlledSystemBase<Scalar, XDimisions, UDimisions> *systemDynamicsPtr() { return systemDynamicsPtr_; }

    int run(const Scalar initTime, const Vector<Scalar, XDimisions> &initState, const Scalar finalTime, ControllerBase<Scalar, XDimisions, UDimisions> *controller,
            std::array<Scalar, ArrayLen> &timeTrajectory, std::array<Vector<Scalar, XDimisions>, ArrayLen> &stateTrajectory, std::array<Vector<Scalar, UDimisions>, ArrayLen> &inputTrajectory) override
    {
        assert(finalTime > initTime);
        assert(controller != nullptr);

        // set controller
        systemDynamicsPtr_->setController(controller);

        Observer<Scalar, XDimisions> observer(ArrayLen, stateTrajectory.data(), timeTrajectory.data()); // concatenate trajectory
        // integrate controlled system
        RK45Intergraor_.integrateConst(*systemDynamicsPtr_, observer, initState, initTime, finalTime, this->settings().timeStep);

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
    ControlledSystemBase<Scalar, XDimisions, UDimisions> *systemDynamicsPtr_{nullptr};

    RungeKuttaDormandPrince5<Scalar, XDimisions> RK45Intergraor_;
};
