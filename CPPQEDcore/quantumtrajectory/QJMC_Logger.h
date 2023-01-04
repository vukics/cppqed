// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Archive.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>

#include <iosfwd>
#include <list>


namespace quantumtrajectory::qjmc {


class Logger;

using LogLevel=std::bitset<3>;


struct EnsembleLogger
{
  typedef std::list<Logger> LoggerList;

  /// Called by Ensemble<QuantumJumpMonteCarlo>::logOnEnd, it calculates a temporal histogram of total jumps
  /** \todo The different kinds of jumps should be collected into different histograms */
  static std::ostream& stream(std::ostream&, const LoggerList&, size_t nBins, size_t nJumpsPerBin);

  template <typename SingleTrajectory>
  auto& operator()(const std::vector<SingleTrajectory>& trajs, std::ostream& os) const
  {
    LoggerList loggerList;
    for (auto& traj : trajs) loggerList.push_back(traj.getLogger());
  
    return stream(os,loggerList,nBins_,nJumpsPerBin_);
  }
  
  const size_t nBins_, nJumpsPerBin_;
  
};



/// Essentially an aggregate of data fields for logging during a QuantumJumpMonteCarlo run.
/**
 * Different log levels are to be supplied as the first parameter to Logger::Logger, and mean the following:
 * 
 * logLevel | effect
 * -------- | ------
 * <=1      | No log output *during* the trajectory.
 * >0       | Summary log information at the end of the trajectory stream
 * >1       | Reporting jumps on `std::cout` also during trajectory evolution.
 * >2       | Reporting `dpLimit` overshoots and the resulting stepsize decrease also during trajectory evolution.
 * >3       | Reporting number of failed ODE steps in the given step of the trajectory evolution.
 * 
 */
class Logger
{
public:
  typedef std::list<std::pair<double,size_t> > jumpTrajectory; ///< Stores <time instant, lindbladNo> pairs, that is, the complete stochastic MCWF trajectory (the rest is deterministic)

  /// Straightforward constructor
  Logger(int logLevel, ///< governs amount of log output during trajectory stream
              size_t nLindblads ///< number of Lindblad operators
             );

  void step(); ///< registers an MCWF step

  void stepBack(std::ostream&, double dp, double    dtDid, double newDtTry, double t, bool logControl); ///< registers a step-back upon \link qjmc::Pars::overshootTolerance tolerance-overshoot\endlink
  void overshot(std::ostream&, double dp, double oldDtTry, double newDtTry          , bool logControl); ///< registers a \link qjmc::Pars::dpLimit dpLimit\endlink overshoot

  void processNorm(double norm); ///< bookkeeps maximal deviation of the stochastic state vector from norm 1

  void jumpOccured(std::ostream&, double t, size_t lindbladNo); ///< registers a jump at time `t` with its identifying ordinal

  std::ostream& onEnd(std::ostream&) const; ///< streams summary log information at the end (called by QuantumJumpMonteCarlo::logOnEnd_v)
  
  const jumpTrajectory& getTrajectory() const {return traj_;}
  
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & nMCWF_steps_ & nOvershot_ & nToleranceOvershot_
                                                      & dpMaxOvershoot_ & dpToleranceMaxOvershoot_ & normMaxDeviation_
                                                      & traj_;}

  friend struct EnsembleLogger;
  
  const int logLevel_;
  const size_t nLindblads_;

  size_t nMCWF_steps_, nOvershot_, nToleranceOvershot_;
  double dpMaxOvershoot_, dpToleranceMaxOvershoot_, normMaxDeviation_;

  jumpTrajectory traj_;
  
};


} // quantumtrajectory::qjmc

