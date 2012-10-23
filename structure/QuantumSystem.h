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
  typedef boost::shared_ptr<const QuantumSystem> Ptr;

  typedef DimensionsBookkeeper<RANK> Base;

  typedef typename Base::Dimensions Dimensions;

  using Base::getDimensions;

  explicit QuantumSystem(const Dimensions& dimensions) : Base(dimensions) {}

  virtual ~QuantumSystem() {}

  double highestFrequency (                ) const {return  highestFrequency_v(  );}
  void   displayParameters(std::ostream& os) const {return displayParameters_v(os);}

private:
  virtual double  highestFrequency_v(             ) const = 0;
  virtual void   displayParameters_v(std::ostream&) const = 0;

};


} // structure

#endif // STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
