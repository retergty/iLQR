#include "iLQR.hpp"

#include "DefaultInitializer.hpp"

int main() {
  DefaultInitializer<double, 3, 2> init;
  DDPSettings<double> ddp_setting;
  using Descriptor = iLQRDescriptor<double, Dimensions<3, 2, 10>>;
  iLQR<Descriptor> ilqr(ddp_setting, nullptr, &init);
  decltype(ilqr)::StateVector_t init_state;
  ilqr.run(0, init_state);
  return 0;
}