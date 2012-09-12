// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

#include "ParsStochasticTrajectory.h"


namespace quantumtrajectory {


struct ParsMCWF_Trajectory : public trajectory::ParsStochasticTrajectory {
  
  double &dpLimit, &overshootTolerance;

  unsigned &svdc;
  bool &firstSVDisplay;
  int &svdPrecision;

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  bool &binarySVFile;
#endif // DO_NOT_USE_BOOST_SERIALIZATION

  std::string &initFile;

  size_t &basisDim;
  std::string &basisFile;

  int &logLevel;

  ParsMCWF_Trajectory(parameters::ParameterTable& p, const std::string& mod="");

};


} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_PARSMCWF_TRAJECTORY_H_INCLUDED
