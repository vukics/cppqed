#include "ParsEvolution.h"

#include "Evolution.h"

#include "Pars.h"


ParsEvolution::ParsEvolution(parameters::ParameterTable& p, const std::string& mod) 
  : ParsMCWF_Trajectory(p,mod),
    evol(p.addTitle("Evolution",mod).addMod("evol",mod,"Evolution mode (single, ensemble, master, convergence)",EM_SINGLE)),
    negativity(p.addMod("negativity",mod,"Calculates negativity in ensemble & master",false))
{}


