// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_PARSSTOCHASTICTRAJECTORY_H_INCLUDED


#include "StochasticTrajectoryFwd.h"

#include "ParsTrajectory.h"


namespace trajectory {


/// Aggregate of parameters pertaining to stochastic simulations
/** \copydetails ParsRun */
struct ParsStochastic : ParsEvolved
{

  unsigned long &seed; ///< random-number generator seed
  
  /// whether the noise should be on or off
  /**
   * (if it makes sense to turn it off at all for a concrete Stochastic
   * – e.g. for a \link quantumtrajectory::MCWF_Trajectory Monte Carlo wave-function trajectory\endlink, turning off the noise means simply to disable quantum jumps)
   */
  bool &noise;
  
  size_t &nTraj; ///< number of trajectories in case of ensemble averaging

  ParsStochastic(parameters::ParameterTable&, const std::string& mod="");
      
};


} // trajectory


#endif // CPPQEDCORE_UTILS_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
