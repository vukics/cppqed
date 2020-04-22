// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED

#include "Archive.h"

#include "core_config.h"

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION

#include <iosfwd>
#include <list>


namespace quantumtrajectory {

  
namespace mcwf {

class Logger;

} // mcwf

/// Auxiliary tools to EnsembleMCWF
namespace ensemble {

using namespace mcwf;
  
typedef std::list<Logger> LoggerList;

/// Called by EnsembleMCWF::logOnEnd, it displays a time histogram of total jumps
/** \todo The different kinds of jumps should be collected into different histograms */
std::ostream& displayLog(std::ostream&, const LoggerList&, size_t nBins, size_t nJumpsPerBin);
  
} // ensemble


namespace mcwf {

/// Essentially an aggregate of data fields for logging during a MCWF_Trajectory run.
/**
 * Different log levels are to be supplied as the first parameter to Logger::Logger, and mean the following:
 * 
 * logLevel | effect
 * -------- | ------
 * <=1      | No log output *during* the trajectory.
 * >0       | Summary log information at the end of the trajectory display
 * >1       | Reporting jumps on `std::cout` also during trajectory evolution.
 * >2       | Reporting `dpLimit` overshoots and the resulting stepsize decrease also during trajectory evolution.
 * >3       | Reporting number of failed ODE steps in the given step of the trajectory evolution.
 * 
 */
class Logger
{
public:
  typedef std::list<std::pair<double,size_t> > MCWF_Trajectory; ///< Stores <time instant, lindbladNo> pairs, that is, the complete stochastic MCWF trajectory (the rest is deterministic)

  /// Straightforward constructor
  Logger(int logLevel, ///< governs amount of log output during trajectory display
              size_t nLindblads ///< number of Lindblad operators
             );

  void step(); ///< registers an MCWF step

  void stepBack(std::ostream&, double dp, double    dtDid, double newDtTry, double t, bool logControl); ///< registers a step-back upon \link mcwf::Pars::overshootTolerance tolerance-overshoot\endlink
  void overshot(std::ostream&, double dp, double oldDtTry, double newDtTry          , bool logControl); ///< registers a \link mcwf::Pars::dpLimit dpLimit\endlink overshoot

  void processNorm(double norm); ///< bookkeeps maximal deviation of the stochastic state vector from norm 1

  void jumpOccured(std::ostream&, double t, size_t lindbladNo); ///< registers a jump at time `t` with its identifying ordinal

  std::ostream& onEnd(std::ostream&) const; ///< displays summary log information at the end (called by MCWF_Trajectory::logOnEnd_v)
  
  const MCWF_Trajectory& getTrajectory() const {return traj_;}
  
private:
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & nMCWF_steps_ & nOvershot_ & nToleranceOvershot_
                                                      & dpMaxOvershoot_ & dpToleranceMaxOvershoot_ & normMaxDeviation_
                                                      & traj_;}
#endif // DO_NOT_USE_BOOST_SERIALIZATION

  friend std::ostream& ensemble::displayLog(std::ostream&, const ensemble::LoggerList&, size_t, size_t);
  
  const int logLevel_;
  const size_t nLindblads_;

  size_t nMCWF_steps_, nOvershot_, nToleranceOvershot_;
  double dpMaxOvershoot_, dpToleranceMaxOvershoot_, normMaxDeviation_;

  MCWF_Trajectory traj_;
  
};


} } // quantumtrajectory::mcwf


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED
