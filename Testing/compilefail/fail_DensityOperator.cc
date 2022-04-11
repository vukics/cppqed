// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "DensityOperator.h"

auto f(quantumdata::DensityOperator<5> rho)
{
  return rho(1,2,3,4,1)(2,1,5,3);
}
