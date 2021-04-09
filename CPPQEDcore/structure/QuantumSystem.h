// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_QUANTUMSYSTEM_H_INCLUDED

#include "DimensionsBookkeeper.h"

#include <memory>


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
private:
  typedef DimensionsBookkeeper<RANK> Base;

public:
  typedef typename Base::Dimensions Dimensions;

  explicit QuantumSystem(const Dimensions& dimensions) : Base(dimensions) {} ///< Construction from a set of Dimensions

  virtual ~QuantumSystem() {}

  double        highestFrequency (                ) const {return  highestFrequency_v(  );} ///< The fastest timescale of the system for ODE stepping
  std::ostream& streamParameters(std::ostream& os) const {return streamParameters_v(os);} ///< Communicating system parameters towards the user

private:
  virtual double         highestFrequency_v(            ) const = 0;
  virtual std::ostream& streamParameters_v(std::ostream&) const = 0;

};


template <int RANK>
using QuantumSystemPtr=std::shared_ptr<const QuantumSystem<RANK>>;
  



} // structure

#endif // CPPQEDCORE_STRUCTURE_QUANTUMSYSTEM_H_INCLUDED
