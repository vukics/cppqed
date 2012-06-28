// -*- C++ -*-
#ifndef _PARS_MCWF_TRAJECTORY_H
#define _PARS_MCWF_TRAJECTORY_H

#include "MCWF_TrajectoryFwd.h"

#include "ParsStochasticTrajectory.h"

#include "EvolvedGSL.h"

namespace quantumtrajectory {


struct ParsMCWF_Trajectory : public trajectory::ParsStochasticTrajectory {
  
  double &dpLimit, &overshootTolerance;

  unsigned &svdc;
  bool &firstSVDisplay;
  int &svdPrecision;

#ifdef USE_BOOST_SERIALIZATION
  bool &binarySVFile;
#endif // USE_BOOST_SERIALIZATION

  std::string &initFile;

  size_t &basisDim;
  std::string &basisFile;

  int& logLevel;

  evolved::SteppingFunction& sf;

  ParsMCWF_Trajectory(parameters::ParameterTable& p, const std::string& mod="");

};


} // quantumtrajectory

#endif // _PARS_MCWF_TRAJECTORY_H
