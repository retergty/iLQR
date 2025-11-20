#include <iostream>
#include "iLQR.hpp"

int main()
{
    iLQR<double, 3, 2, 10> ilqr;
    decltype(ilqr)::StateVector_t init_state;
    ilqr.run(0, init_state, 10);
    return 0;
}