// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ExpectationValues.h"

#include "FormDouble.h"


std::ostream& structure::stream(const EV_Array& eva, std::ostream& os, int precision)
{
  os<<'\t';
  {
    const FormDouble fd{precision};
    for (double ev : eva) os<<fd(ev);
  }
  return os;
}


void structure::calculateVariance(EV_Array& eva)
{
  eva[1]-=sqr(eva[0]);
}
