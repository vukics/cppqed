// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED
#define QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED

#include <iosfwd>
#include <list>


namespace quantumtrajectory {


class MCWF_TrajectoryLogger
{
public:
  typedef std::list<std::pair<double,size_t> > MCWF_Trajectory;
  // Stores <time instant, jumpNo> pairs

  MCWF_TrajectoryLogger(unsigned logLevel, bool isHamiltonian, std::ostream& os);

  ~MCWF_TrajectoryLogger();

  void step() const;

  void stepBack(double dp, double    dtDid, double newDtTry, double t, bool logControl) const;
  void overshot(double dp, double oldDtTry, double newDtTry          , bool logControl) const;

  void processNorm(double norm) const;

  void jumpOccured(double t, size_t jumpNo) const;

  void logFailedSteps(size_t) const;

  void hamiltonianCalled() const;

private:
  const int logLevel_;
  const bool isHamiltonian_;
  std::ostream& os_;

  mutable size_t nSteps_, nOvershot_, nToleranceOvershot_, nFailedSteps_, nHamiltonianCalls_;
  mutable double dpMaxOvershoot_, dpToleranceMaxOvershoot_, normMaxDeviation_;

  mutable MCWF_Trajectory traj_;
  
};


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MCWF_TRAJECTORYLOGGER_H_INCLUDED
