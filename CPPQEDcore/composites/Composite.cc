// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "Composite.h"

#include "SubSystem.h"

namespace composite {

bool compareFreesFrequency(const SubSystemFree& ssf1, const SubSystemFree& ssf2)
{
  return ssf1.get()->highestFrequency() < ssf2.get()->highestFrequency();
}

} // composite
