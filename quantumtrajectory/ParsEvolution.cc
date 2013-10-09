#include "ParsEvolution.h"

#include "Evolution_.h"

#include "Pars.tcc"


ParsEvolution::ParsEvolution(parameters::ParameterTable& p, const std::string& mod) 
  : ParsRun(p,mod),
    ParsMCWF(p,mod),
    evol(p.addTitle("Evolution",mod).addMod("evol",mod,"Evolution mode (single, ensemble, master, master_fast)",EM_SINGLE)),
    negativity(p.addMod("negativity",mod,"Calculates negativity in ensemble & master",false)),
    timeAverage(p.addMod("timeAverage",mod,"Calculates time averages in MCWF trajectory",false)),
    relaxationTime(p.addMod("relaxationTime",mod,"Relaxation time for time averaging",0.))
{}


