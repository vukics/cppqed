// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_NONORTHOGONALDENSITYOPERATOR_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_NONORTHOGONALDENSITYOPERATOR_H_INCLUDED

#include "NonOrthogonalDensityOperatorFwd.h"

#include "DensityOperator.h"


namespace quantumdata {

template<int RANK, typename TRAFO>
class NonOrthogonalDensityOperator/* : public DensityOperator<RANK>, private linalg::VectorSpace<NonOrthogonalDensityOperator<RANK,TRAFO> > {

  typedef StateVectorBase<RANK> Base;
  typedef Base::Impl Impl;
  typedef Base::Idx  Idx ;

  NonOrthogonalStateVector(const Impl& array) : Base(array), _arrayPulled(array.shape()), _PulledCompatible(false) {}
  // Not pulled in the first place

  Impl& operator()() {_PulledCompatible=false; return const_cast<Impl&>(static_cast<const Base*>(this)->operator());}

  // matrix-like indexing
  const dcomp matrix(const Idx& i1, const Idx& i2);

  renorm

  CLDO MakeLDO() const; // For this it has to pull

private:
  Impl _arrayPulled;
  bool _PulledCompatible;

  strategyForPulling _p; 

};

*/;

} // quantumdata


#endif // CPPQEDCORE_QUANTUMDATA_NONORTHOGONALDENSITYOPERATOR_H_INCLUDED
