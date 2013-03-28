// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED
#define QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

#include "Archive.h"

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION

#include <iosfwd>
#include <list>


namespace quantumtrajectory {


namespace ensemblemcwf {

typedef std::list<const MCWF_Logger*> LoggerList;

std::ostream& displayLog(std::ostream&, const LoggerList&);
  
} // ensemblemcwf


class MCWF_Logger
{
public:
  typedef std::list<std::pair<double,size_t> > MCWF_Trajectory;
  // Stores <time instant, jumpNo> pairs

  MCWF_Logger(int logLevel, bool isHamiltonian, size_t nJumps);

  void step() const;

  void stepBack(double dp, double    dtDid, double newDtTry, double t, bool logControl) const;
  void overshot(double dp, double oldDtTry, double newDtTry          , bool logControl) const;

  void processNorm(double norm) const;

  void jumpOccured(double t, size_t jumpNo) const;

  void logFailedSteps(size_t) const;

  void hamiltonianCalled() const;

  std::ostream& onEnd(std::ostream&) const;
  
private:
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & nSteps_ & nOvershot_ & nToleranceOvershot_ & nFailedSteps_ & nHamiltonianCalls_
                                                      & dpMaxOvershoot_ & dpToleranceMaxOvershoot_ & normMaxDeviation_
                                                      & traj_;}
#endif // DO_NOT_USE_BOOST_SERIALIZATION

  friend std::ostream& ensemblemcwf::displayLog(std::ostream&, const ensemblemcwf::LoggerList&);
  
  const int logLevel_;
  const bool isHamiltonian_;
  const size_t nJumps_;

  mutable size_t nSteps_, nOvershot_, nToleranceOvershot_, nFailedSteps_, nHamiltonianCalls_;
  mutable double dpMaxOvershoot_, dpToleranceMaxOvershoot_, normMaxDeviation_;

  mutable MCWF_Trajectory traj_;
  
};


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED
