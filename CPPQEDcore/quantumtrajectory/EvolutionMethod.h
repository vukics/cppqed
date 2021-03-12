// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Highest level driver functions for quantum trajectories}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_METHOD_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_METHOD_H_INCLUDED

#include <iosfwd>

/// Auxiliary tools for the evolve functions
namespace evolution {

/// Method of evolution for a quantum system
enum Method {
  SINGLE, ///< single \link quantumtrajectory::MCWF_Trajectory MCWF trajectory\endlink
  ENSEMBLE, ///< \link quantumtrajectory::EnsembleMCWF ensemble\endlink of MCWF trajectories
  MASTER ///< Master equation with \link quantumtrajectory::master::Base normal iteration\endlink
};

std::ostream& operator<<(std::ostream&, Method ); ///< output streaming for Method
std::istream& operator>>(std::istream&, Method&); ///< input streaming for Method

} // evolution


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_METHOD_H_INCLUDED
