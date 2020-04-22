// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED

#include "ParsStochasticTrajectory.h"


namespace quantumtrajectory {


/// Auxiliary tools to MCWF_Trajectory  
namespace mcwf {


/// Aggregate of parameters pertaining to \link MCWF_Trajectory MCWF\endlink simulations
/** \copydetails trajectory::ParsRun */
struct Pars : public trajectory::ParsStochastic {
  
  double
    &dpLimit, ///< the parameter \f$\delta p_\text{limit}\f$ (cf. 2.b.ii \link MCWF_Trajectory here\endlink)
    &overshootTolerance; ///< the parameter \f$\delta p_\text{limit}'/\delta p_\text{limit}\f$ (cf. 2.b.ii \link MCWF_Trajectory here\endlink)

  size_t
    &nBins, ///< governs how many bins should be used for the histogram of jumps created by ensemble::displayLog (a zero value means a heuristic automatic determination)
    &nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  Pars(parameters::ParameterTable& p, const std::string& mod="");

};


} } // quantumtrajectory::mcwf

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
