// -*- C++ -*-
#ifndef MCWF_TRAJECTORY_LOGGER_H
#define MCWF_TRAJECTORY_LOGGER_H

#include <iosfwd>
#include <list>


namespace quantumtrajectory {


class MCWF_TrajectoryLogger
{
public:
  typedef std::list<std::pair<double,size_t> > MCWF_Trajectory;
  // Stores <time instant, jumpNo> pairs

  MCWF_TrajectoryLogger(unsigned logLevel, std::ostream& os);

  ~MCWF_TrajectoryLogger();

  void step() const;

  void overshot(double dp, double oldDtTry, double newDtTry) const;

  void processNorm(double norm) const;

  void jumpOccured(double t, size_t jumpNo) const;

private:
  const int logLevel_;
  std::ostream& os_;

  mutable size_t nSteps_, nOvershot_;
  mutable double dpMaxOvershoot_, normMaxDeviation_;

  mutable MCWF_Trajectory traj_;
  
};


} // quantumtrajectory


#endif // MCWF_TRAJECTORY_LOGGER_H
