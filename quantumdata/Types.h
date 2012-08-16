// -*- C++ -*-
#ifndef QUANTUMDATA_TYPES_H_INCLUDED
#define QUANTUMDATA_TYPES_H_INCLUDED

#include "TypesFwd.h"

#include "LazyDensityOperatorFwd.h"

#include "BlitzArray.h"


namespace structure {

struct StaticTag {};

const StaticTag theStaticOne=StaticTag();

} //structure


namespace quantumdata {


namespace details {

struct Empty {};

} // details


template<int RANK, typename B>
struct Types : B 
// it's basically just a metafunction
{
  typedef TTD_CARRAY(  RANK)     StateVectorLow;
  typedef TTD_CARRAY(2*RANK) DensityOperatorLow;
};


} // quantumdata

#endif // QUANTUMDATA_TYPES_H_INCLUDED
