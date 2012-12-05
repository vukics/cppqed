// -*- C++ -*-
#ifndef UTILS_INCLUDE_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
#define UTILS_INCLUDE_PARSSTOCHASTICTRAJECTORY_H_INCLUDED


#include "StochasticTrajectoryFwd.h"

#include "ParsFwd.h"

#include "ParsTrajectory.h"


namespace trajectory {


struct ParsStochastic : Pars {

  unsigned long &seed;
  bool &noise;
  size_t &nTraj;

  ParsStochastic(parameters::ParameterTable&, const std::string& mod="");
      
};


} // trajectory


#endif // UTILS_INCLUDE_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
