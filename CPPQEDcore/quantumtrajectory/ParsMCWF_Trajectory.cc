#include "ParsMCWF_Trajectory.h"

#include "Pars.tcc"


quantumtrajectory::mcwf::Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : ParsStochastic(p,mod),
    dpLimit(p.addTitle("MCWF_Trajectory",mod).addMod("dpLimit",mod,"MCWFS stepper total jump probability limit",0.01)),
    overshootTolerance(p.addMod("overshootTolerance",mod,"Jump probability overshoot tolerance factor",10.)),
    logLevel(p.addMod("logLevel",mod,"logging level",0))
{}

