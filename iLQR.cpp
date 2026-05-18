#include "iLQR.hpp"
#include <iostream>

int main() {
  DefaultInitializer<double, 3, 2> init;
  DDPSettings<double> ddp_setting;
  iLQR<double, 3, 2, 10> ilqr(ddp_setting, nullptr, &init);
  decltype(ilqr)::StateVector_t init_state;
  ilqr.run(0, init_state);
  return 0;
}