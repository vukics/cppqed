// -*- C++ -*-

#include "Composite.h"

#include "SubSystem.h"

namespace composite {

bool compareFreesFrequency(const structure::SubSystemFree& ssf1, const structure::SubSystemFree& ssf2)
{
  return ssf1.get()->highestFrequency() < ssf2.get()->highestFrequency();
}

} // composite
