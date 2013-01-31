// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

#include "ParsStochasticTrajectory.h"


namespace quantumtrajectory {


struct ParsMCWF_Trajectory : public trajectory::ParsStochastic {
  
  double &dpLimit, &overshootTolerance;

  unsigned &svdc;
  bool &firstSVDisplay;
  int &svdPrecision;

  std::string &initFile;

  size_t &basisDim;
  std::string &basisFile;

  int &logLevel;

  ParsMCWF_Trajectory(parameters::ParameterTable& p, const std::string& mod="");

};


} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
