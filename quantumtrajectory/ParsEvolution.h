// -*- C++ -*-
#ifndef _PARS_EVOLUTION_MODE
#define _PARS_EVOLUTION_MODE

#include "Evolution_Fwd.h"

#include "ParsMCWF_Trajectory.h"

using namespace quantumtrajectory;


struct ParsEvolution : public ParsMCWF_Trajectory {

  EvolutionMode &evol;
  bool &negativity;

  ParsEvolution(parameters::ParameterTable& p, const std::string& mod="");

};


#endif // _PARS_EVOLUTION_MODE
