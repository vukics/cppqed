// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines display_densityoperator::_}
// -*- C++ -*-
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED

#include "DensityOperator.h"

#include "Averaged.h"


namespace quantumtrajectory {

/// Contains display_densityoperator::_
namespace display_densityoperator {


/// Wraps common functionality of Master & EnsembleMCWF concerning display of quantum averages on the basis of density operators
/**
 * This comprises
 * - keeping a structure::Averaged instant and calling structure::Averaged::display
 * - performing \link quantumdata::negPT negativity calculation\endlink if needed
 * - extending key with negativity when needed
 * 
 * \tparamRANK
 * \tparam V has the same function as the template parameter `V` in quantumdata::negPT
 * 
 */
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

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED
