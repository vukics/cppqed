// -*- C++ -*-
/// \briefFileDefault
#ifndef QUANTUMDATA_TYPES_H_INCLUDED
#define QUANTUMDATA_TYPES_H_INCLUDED

#include "TypesFwd.h"

#include "LazyDensityOperatorFwd.h"

#include "BlitzArray.h"


namespace quantumdata {


namespace details {

struct Empty {};

} // details


/// Basically only a metafunction defining types for higher-level constructs of arity `RANK`
/**
 * \tparamRANK
 * \tparam B An optional base class for base-class chaining. An empty class (quantumdata::details::Empty) by default
 */
template<int RANK, typename B>
struct Types : B 
// it's basically just a metafunction
{
  typedef TTD_CARRAY(  RANK)     StateVectorLow;
  typedef TTD_CARRAY(2*RANK) DensityOperatorLow;
};


} // quantumdata

#endif // QUANTUMDATA_TYPES_H_INCLUDED
