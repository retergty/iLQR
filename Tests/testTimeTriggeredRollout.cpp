/**
 * @file testTimeTriggeredRollout.cpp
 * @brief TimeTriggeredRollout 兼容性测试：匹配当前项目中的 rollout API。
 */
#include <gtest/gtest.h>

#include <array>

#include "Controller/LinearController.hpp"
#include "Dynamics/LinearSystemDynamics.hpp"
#include "Rollout/TimeTriggeredRollout.hpp"

namespace {
using Scalar = double;
constexpr int XDim = 2;
constexpr int UDim = 1;
constexpr size_t ControllerLen = 2;
constexpr size_t NumIntervals = 1000;
constexpr size_t SampleCount = NumIntervals + 1;

using StateVector = Vector<Scalar, XDim>;
using InputVector = Vector<Scalar, UDim>;
using StateMatrix = Matrix<Scalar, XDim, XDim>;
using InputMatrix = Matrix<Scalar, XDim, UDim>;
using Controller = LinearController<Scalar, XDim, UDim, ControllerLen>;
using System = LinearSystemDynamics<Scalar, XDim, UDim>;
using Rollout = TimeTriggeredRollout<Scalar, XDim, UDim>;
}  // namespace

TEST(TimeTriggeredRolloutCompatibilityTest,
     ProducesConsistentTrajectoryBuffers) {
  constexpr Scalar initTime = 0.0;
  constexpr Scalar finalTime = 10.0;
  constexpr Scalar timeStep = 0.01;

  const StateMatrix A{{-2.0, -1.0}, {1.0, 0.0}};
  const InputMatrix B{1.0, 0.0};
  System systemDynamics(A, B);

  const std::array<Scalar, ControllerLen> controllerTime = {initTime,
                                                            finalTime};
  const InputVector unitInput = InputVector::Ones();
  const std::array<InputVector, ControllerLen> controllerBias = {unitInput,
                                                                 unitInput};
  const Matrix<Scalar, UDim, XDim> zeroGain =
      Matrix<Scalar, UDim, XDim>::Zero();
  const std::array<Matrix<Scalar, UDim, XDim>, ControllerLen> controllerGain = {
      zeroGain, zeroGain};
  Controller controller(controllerTime, controllerBias, controllerGain);

  Rollout rollout(&systemDynamics, timeStep);

  std::array<Scalar, SampleCount> timeTrajectory{};
  std::array<StateVector, SampleCount> stateTrajectory{};
  std::array<InputVector, SampleCount> inputTrajectory{};
  RolloutTrajectoryPointer<Scalar, XDim, UDim> trajectory(
      timeTrajectory.data(), stateTrajectory.data(), inputTrajectory.data(),
      static_cast<int>(SampleCount));

  const StateVector initState = StateVector::Zero();
  const int count =
      rollout.run(initTime, initState, finalTime, &controller, trajectory);

  ASSERT_EQ(count, static_cast<int>(SampleCount));
  EXPECT_NEAR(timeTrajectory.front(), initTime, 1e-12);
  EXPECT_NEAR(timeTrajectory[count - 1], finalTime, 1e-12);
  EXPECT_TRUE(stateTrajectory.front().isApprox(initState, 1e-12));

  for (int i = 0; i < count; ++i) {
    EXPECT_NEAR(inputTrajectory[i](0), 1.0, 1e-12) << "i = " << i;
  }

  EXPECT_NEAR(stateTrajectory[count - 1](0), 0.0, 1e-3);
  EXPECT_NEAR(stateTrajectory[count - 1](1), 1.0, 1e-3);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
