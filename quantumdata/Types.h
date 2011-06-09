// -*- C++ -*-
#ifndef _STRUCTURE_TYPES_H
#define _STRUCTURE_TYPES_H

#include "TypesFwd.h"

#include "LazyDensityOperatorFwd.h"

#include "BlitzArray.h"


namespace structure {

struct StaticTag {};

const StaticTag theStaticOne=StaticTag();


struct LiouvilleanAveragedCommon
{
  typedef TTD_DARRAY(1) DArray1D;
  static const DArray1D defaultArray;
};


} // structure



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

#endif // _STRUCTURE_TYPES_H
