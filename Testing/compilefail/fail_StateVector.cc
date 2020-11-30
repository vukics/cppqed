// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "StateVector.h"

auto f(quantumdata::StateVector<5> psi)
{
  return psi(1,2,3,4,1,0);
}
