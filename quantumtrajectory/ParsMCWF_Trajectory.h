// -*- C++ -*-
#ifndef _PARS_MCWF_TRAJECTORY_H
#define _PARS_MCWF_TRAJECTORY_H

#include "MCWF_TrajectoryFwd.h"

#include "ParsStochasticTrajectory.h"

namespace quantumtrajectory {


struct ParsMCWF_Trajectory : public trajectory::ParsStochasticTrajectory {
  
  double &dpLimit;
  unsigned &svdc;

  std::string &initFile;

  size_t &basisDim;
  std::string &basisFile;

  bool& doLog;

  ParsMCWF_Trajectory(parameters::ParameterTable& p, const std::string& mod="");

};


} // quantumtrajectory

#endif // _PARS_MCWF_TRAJECTORY_H
