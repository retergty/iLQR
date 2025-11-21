#include <iostream>
#include "iLQR.hpp"

int main()
{
    iLQR<double, 3, 2, 10> ilqr(nullptr, 0.1);
    decltype(ilqr)::StateVector_t init_state;
    ilqr.run(0, init_state);
    return 0;
}