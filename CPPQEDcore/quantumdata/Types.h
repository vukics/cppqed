// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED

#include "BlitzArray.h"


namespace quantumdata {


namespace details { struct EmptyBase {}; }

/// Basically only a metafunction defining types for higher-level constructs of arity `RANK`
/**
 * \tparamRANK
 * \tparam B An optional base class for base-class chaining. An empty class by default
 */
template<int RANK, typename B=details::EmptyBase>
struct Types : B 
// it's basically just a metafunction
{
  typedef CArray<  RANK>     StateVectorLow;
  typedef CArray<2*RANK> DensityOperatorLow;
};


} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_TYPES_H_INCLUDED
