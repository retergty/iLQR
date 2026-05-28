#include "iLQR/iLQR.hpp"

#include "Dynamics/LinearSystemDynamics.hpp"
#include "Initialization/DefaultInitializer.hpp"

int main() {
  DefaultInitializer<double, 3, 2> init;
  DDPSettings<double> ddp_setting;
  using Descriptor =
      iLQRDescriptor<double,
                     TranscriptionConfig<Dimensions<3, 2>, Horizon<10>>>;
  using Solver = iLQR<Descriptor>;
  Matrix<double, 3, 3> A = Matrix<double, 3, 3>::Zero();
  Matrix<double, 3, 2> B = Matrix<double, 3, 2>::Zero();
  LinearSystemDynamics<double, 3, 2> dynamics(A, B);
  Solver::OptimalControlProblem_t problem;
  problem.dynamicsPtr = &dynamics;
  Solver ilqr(ddp_setting, problem, &init);
  decltype(ilqr)::StateVector_t init_state;
  init_state.setZero();
  ilqr.run(0, init_state);
  return 0;
}