// -*- C++ -*-
#ifndef UTILS_INCLUDE_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
#define UTILS_INCLUDE_PARSSTOCHASTICTRAJECTORY_H_INCLUDED


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


#endif // UTILS_INCLUDE_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
