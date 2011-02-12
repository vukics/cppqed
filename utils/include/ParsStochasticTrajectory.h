// -*- C++ -*-
#ifndef _PARS_STOCHASTIC_TRAJECTORY_H
#define _PARS_STOCHASTIC_TRAJECTORY_H


#include "StochasticTrajectoryFwd.h"

#include "ParsFwd.h"

#include "ParsTrajectory.h"


namespace trajectory {


struct ParsStochasticTrajectory : ParsTrajectory {

  unsigned long &seed;
  bool &noise;
  size_t &nTraj;

  ParsStochasticTrajectory(parameters::ParameterTable&, const std::string& mod="");
      
};


} // trajectory


#endif // _PARS_STOCHASTIC_TRAJECTORY_H
