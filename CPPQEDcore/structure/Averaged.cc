// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Averaged.h"

#include "FormDouble.h"


std::ostream& stream(const EV_Array& eva, std::ostream& os, int precision)
{
  os<<'\t';
  {
    const formdouble::FormDouble fd{precision};
    for (double ev : eva) os<<fd(ev);
  }
  return os;
}
