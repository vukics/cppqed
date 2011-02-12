// -*- C++ -*-
#ifndef _STRUCTURE_QUANTUM_SYSTEM_H
#define _STRUCTURE_QUANTUM_SYSTEM_H

#include "QuantumSystemFwd.h"

#include "DimensionsBookkeeper.h"


namespace structure {

// The QuantumSystem class is the abstract interface every system has
// to present towards the drivers MCWF_Trajectory, Ensemble, Master


template<int RANK>
class QuantumSystem : public DimensionsBookkeeper<RANK>
{
public:
  typedef DimensionsBookkeeper<RANK> Base;

  typedef typename Base::Dimensions Dimensions;

  using Base::getDimensions;

  explicit QuantumSystem(const Dimensions& dimensions) : Base(dimensions) {}

  virtual ~QuantumSystem() {}

  virtual double highestFrequency (             ) const = 0;
  virtual void   displayParameters(std::ostream&) const = 0;

};


} // structure

#endif // _STRUCTURE_QUANTUM_SYSTEM_H
