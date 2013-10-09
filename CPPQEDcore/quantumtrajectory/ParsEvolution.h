// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_PARSEVOLUTION_H_INCLUDED
#define QUANTUMTRAJECTORY_PARSEVOLUTION_H_INCLUDED

#include "Evolution_Fwd.h"

#include "ParsMCWF_Trajectory.h"

using namespace quantumtrajectory;


struct ParsEvolution : public trajectory::ParsRun, public ParsMCWF {

  EvolutionMode &evol;
  bool &negativity, &timeAverage;
  double &relaxationTime;

  ParsEvolution(parameters::ParameterTable& p, const std::string& mod="");

};


#endif // QUANTUMTRAJECTORY_PARSEVOLUTION_H_INCLUDED
