// -*- C++ -*-
#ifndef UTILS_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
#define UTILS_PARSSTOCHASTICTRAJECTORY_H_INCLUDED


#include "StochasticTrajectoryFwd.h"

#include "ParsTrajectory.h"


namespace trajectory {


struct ParsStochastic : ParsEvolved
{

  unsigned long &seed;
  bool &noise;
  size_t &nTraj;

  ParsStochastic(parameters::ParameterTable&, const std::string& mod="");
      
};


} // trajectory


#endif // UTILS_PARSSTOCHASTICTRAJECTORY_H_INCLUDED
