// -*- C++ -*-
#ifndef STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
#define STRUCTURE_QUANTUMSYSTEM_H_INCLUDED

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

#endif // STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
