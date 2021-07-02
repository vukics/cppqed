// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED

#include "BlitzArray.h"


namespace structure {

using Averages = DArray<1>; using Rates = Averages;

} // structure


namespace quantumdata {

template<int RANK>
using StateVectorLow=CArray<RANK>;

template<int RANK>
using DensityOperatorLow=CArray<2*RANK>;
  

} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED
