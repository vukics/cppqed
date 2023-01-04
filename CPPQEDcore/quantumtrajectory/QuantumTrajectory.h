// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"

#include "Trajectory.h"

#include <bitset>


namespace quantumtrajectory {

using EntanglementMeasuresSwitch = std::bitset<3>;


using StreamReturnType=std::tuple<std::ostream&,::structure::EV_Array>;

  
/// Forwards to trajectory::initialTimeStep, with the highest frequency of the system taken as structure::QuantumSystem::highestFrequency
template<size_t RANK, ::structure::quantum_system_dynamics<RANK> QSD>
inline double initialTimeStep(structure::QuantumSystemPtr<RANK> qs)
{
  return cppqedutils::trajectory::initialTimeStep(qs->highestFrequency());
}


template<size_t RANK, ::structure::quantum_system_dynamics<RANK> QSD>
std::ostream& streamCharacteristics(const QSD& qsd, std::ostream& os)
{
  return os<<"System characteristics: "
      <<(castEx(qs) ? "Interaction picture, "   : "")
      <<(castHa(qs) ? "Hamiltonian evolution, " : "")
      <<(castLi(qs) ? "Liouvillian evolution, " : "")
      <<(castAv(qs) ? "calculates Averages."    : "");
}


} // quantumtrajectory

