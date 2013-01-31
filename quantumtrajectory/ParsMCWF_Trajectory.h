// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

#include "ParsStochasticTrajectory.h"


namespace quantumtrajectory {


struct ParsMCWF : public trajectory::ParsStochastic {
  
  double &dpLimit, &overshootTolerance;

  int &logLevel;

  ParsMCWF(parameters::ParameterTable& p, const std::string& mod="");

};


} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
