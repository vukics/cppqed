// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
#define STRUCTURE_QUANTUMSYSTEM_H_INCLUDED

#include "QuantumSystemFwd.h"

#include "DimensionsBookkeeper.h"

#include <boost/shared_ptr.hpp>


namespace structure {

/// The abstract interface every system has to present towards the quantum trajectory drivers quantumtrajectory::MCWF_Trajectory, quantumtrajectory::EnsembleMCWF, and quantumtrajectory::Master
/**
 * Describes an entity that has dimensions in a Hilbert space of arity RANK, it may have some frequencies, and some parameters to communicate towards the user.
 * 
 * \tparamRANK
 * 
 * \note The class features the non-virtual interface idiom, which is ubiquitous in the framework.
 */
template<int RANK>
class QuantumSystem : public DimensionsBookkeeper<RANK>
{
public:
  typedef boost::shared_ptr<const QuantumSystem> Ptr; ///< Many of the basic template classes in the framework act as template metafunctions returning a shared pointer to their own type
  
private:
  typedef DimensionsBookkeeper<RANK> Base;

public:
  typedef typename Base::Dimensions Dimensions;

  using Base::getDimensions;

  explicit QuantumSystem(const Dimensions& dimensions) : Base(dimensions) {} ///< Construction from a set of Dimensions

  virtual ~QuantumSystem() {}

  double        highestFrequency (                ) const {return  highestFrequency_v(  );} ///< The fastest timescale of the system for ODE stepping
  std::ostream& displayParameters(std::ostream& os) const {return displayParameters_v(os);} ///< Communicating system parameters towards the user

private:
  virtual double         highestFrequency_v(             ) const = 0;
  virtual std::ostream& displayParameters_v(std::ostream&) const = 0;

};


} // structure

#endif // STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
