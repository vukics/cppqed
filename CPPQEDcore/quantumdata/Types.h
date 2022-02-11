// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED

#include "BlitzArray.h"


namespace quantumdata {

template<int RANK>
using StateVectorLow=CArray<RANK>;

template<int RANK>
using DensityOperatorLow=CArray<2*RANK>;
  

} // quantumdata


namespace structure {

using ::quantumdata::StateVectorLow, ::quantumdata::DensityOperatorLow;

using Averages = DArray<1>; using Rates = Averages;

} // structure



#endif // CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED
