#include "ParsMCWF_Trajectory.h"

#include "Pars.h"


namespace quantumtrajectory {


ParsMCWF_Trajectory::ParsMCWF_Trajectory(parameters::ParameterTable& p, const std::string& mod)
  : ParsStochasticTrajectory(p,mod),
    dpLimit(p.addTitle("MCWF_Trajectory",mod).addMod("dpLimit",mod,"MCWFS stepper total jump probability limit",0.1)),
    svdc(p.addMod("svdc",mod,"Number of displays between two state-vector Displays",0u)),
    initFile (p.addMod<std::string>("initFile" ,mod,"file containing the initial state vector","")),
    basisDim(p.addMod("basisDim",mod,"number of basis vectors the stochastic wave function is compared against",size_t(0))),
    basisFile(p.addMod<std::string>("basisFile",mod,"file containing the basis vectors","")),
    doLog(p.addMod("doLog",mod,"logging level",0u)),
    sf(p.addMod("steppingFunction",mod,"Stepping function for EvolvedGSL",evolved::SF_RKCK))
{}


} // quantumtrajectory
