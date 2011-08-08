// -*- C++ -*-
#ifndef MCWF_TRAJECTORY_LOGGER_H
#define MCWF_TRAJECTORY_LOGGER_H

#include <iosfwd>


namespace quantumtrajectory {


class MCWF_TrajectoryLogger
{
public:
  MCWF_TrajectoryLogger(unsigned logLevel, std::ostream& os);

  ~MCWF_TrajectoryLogger();

  void step() const;

  void overshot(double dp, double oldDtTry, double newDtTry) const;

  void processNorm(double norm) const;

private:
  const unsigned logLevel_;
  std::ostream& os_;

  mutable size_t nSteps_, nOvershot_;
  mutable double dpMaxOvershoot_, normMaxDeviation_;
  
};


} // quantumtrajectory


#endif // MCWF_TRAJECTORY_LOGGER_H
