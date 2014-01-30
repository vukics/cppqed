// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED
#define QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED

#include "DensityOperator.h"

#include "Averaged.h"


namespace quantumtrajectory {

namespace display_densityoperator {


template<int RANK, typename V>
class _
{
public:
  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  typedef          structure::Averaged<RANK>      Averaged   ;
  typedef typename structure::Averaged<RANK>::Ptr AveragedPtr;

  _(AveragedPtr av, bool negativity) : av_(av), negativity_(negativity) {}

  std::ostream& display   (double t, const DensityOperator&, std::ostream&, int precision) const;
  std::ostream& displayKey(std::ostream&, size_t&) const;

private:
  const AveragedPtr av_ ;

  const bool negativity_;

};


} // display_densityoperator


} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED
