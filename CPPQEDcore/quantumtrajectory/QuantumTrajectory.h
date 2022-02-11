// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED

#include "QuantumSystem.h"
#include "Structure.h"

#include "Trajectory.h"


/// Comprises modules representing trajectory drivers for simulating quantum systems
namespace quantumtrajectory {


using StreamReturnType=std::tuple<std::ostream&,structure::Averages>;

  
/// Forwards to trajectory::initialTimeStep, with the highest frequency of the system taken as structure::QuantumSystem::highestFrequency
template<int RANK>
inline double initialTimeStep(structure::QuantumSystemPtr<RANK> qs)
{
  return cppqedutils::trajectory::initialTimeStep(qs->highestFrequency());
}


} // quantumtrajectory

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
