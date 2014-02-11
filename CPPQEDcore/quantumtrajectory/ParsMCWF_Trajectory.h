/// \briefFileDefault
#ifndef QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

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

  int &logLevel; ///< governs how much logging information is displayed during an MCWF_Trajectory run \see Logger

  Pars(parameters::ParameterTable& p, const std::string& mod="");

};


} } // quantumtrajectory::mcwf

#endif // QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
