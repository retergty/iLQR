#include <iostream>
#include "iLQR.hpp"

int main()
{
    DefaultInitializer<double, 3, 2> init;
    iLQR<double, 3, 2, 10> ilqr(nullptr, &init);
    decltype(ilqr)::StateVector_t init_state;
    ilqr.run(0, init_state);
    return 0;
}