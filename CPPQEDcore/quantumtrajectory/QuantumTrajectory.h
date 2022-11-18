// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED

#include "ExpectationValues.h"
#include "Liouvillian.h"
#include "Hamiltonian.h"
#include "QuantumSystem.h"

#include "Trajectory.h"

#include <bitset>


namespace quantumtrajectory {

using EntanglementMeasuresSwitch = std::bitset<3>;

using StreamReturnType=std::tuple<std::ostream&,structure::Averages>;

  
/// Forwards to trajectory::initialTimeStep, with the highest frequency of the system taken as structure::QuantumSystem::highestFrequency
template<int RANK>
inline double initialTimeStep(structure::QuantumSystemPtr<RANK> qs)
{
  return cppqedutils::trajectory::initialTimeStep(qs->highestFrequency());
}


template<size_t RANK, ::structure::hamiltonian<RANK> HA, ::structure::liouvillian<RANK> LI, ::structure::expectation_values<RANK> EV>
std::ostream& streamCharacteristics(const HA& ha, const LI& li, const EV& ev, std::ostream& os);/*
{
  return os<<"System characteristics: "
      <<(castEx(qs) ? "Interaction picture, "   : "")
      <<(castHa(qs) ? "Hamiltonian evolution, " : "")
      <<(castLi(qs) ? "Liouvillian evolution, " : "")
      <<(castAv(qs) ? "calculates Averages."    : "");
}*/


} // quantumtrajectory

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
