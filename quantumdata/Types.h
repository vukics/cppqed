// -*- C++ -*-
/// \briefFileDefault
#ifndef QUANTUMDATA_TYPES_H_INCLUDED
#define QUANTUMDATA_TYPES_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include "BlitzArray.h"

#include <boost/mpl/empty_base.hpp>


namespace quantumdata {


/// Basically only a metafunction defining types for higher-level constructs of arity `RANK`
/**
 * \tparamRANK
 * \tparam B An optional base class for base-class chaining. An empty class (boost::mpl::empty_base) by default
 */
template<int RANK, typename B=boost::mpl::empty_base>
struct Types : B 
// it's basically just a metafunction
{
  typedef CArray<  RANK>     StateVectorLow;
  typedef CArray<2*RANK> DensityOperatorLow;
};


} // quantumdata

#endif // QUANTUMDATA_TYPES_H_INCLUDED
