#include "ParsStochasticTrajectory.h"

#include "FormDouble.h"
#include "Pars.tcc"

namespace trajectory {


const double ParsEvolved::epsRelDefault=1e-6;
const double ParsEvolved::epsAbsDefault=1e-30;


ParsRun::ParsRun(parameters::ParameterTable& p, const std::string& mod)
  : T(p.addTitle("Trajectory",mod).addMod("T",mod,"Simulated time",1.)),
    dc(p.addMod("dc",mod,"Number of steps between two Displays",10)),
    Dt(p.addMod("Dt",mod,"Timestep between two Displays",.1)),
    NDt(p.addMod("NDt",mod,"Number of steps in Dt mode",0L)),
    ofn(p.addMod<std::string>("o",mod,"Output file name for Trajectory, when empty, cout","")),
    initialFileName(p.addMod<std::string>("initialFileName",mod,"Trajectory initial file name","")),
    precision(p.addMod("precision",mod,"General precision of output",formdouble::Zero(FormDouble::defaultPrecision/2))),
    displayInfo(p.addMod("displayInfo",mod,"Whether to display header for trajectories",true)),
    firstStateDisplay(p.addMod("firstStateDisplay",mod,"Displays state vector at startup",true)),
    sdf(p.addMod("sdf",mod,"State output frequency",0u))
{}


ParsEvolved::ParsEvolved(parameters::ParameterTable& p, const std::string& mod)
  : epsRel(p.addTitle("Evolved",mod).addMod("eps"   ,mod,"ODE stepper relative precision",epsRelDefault)),
    epsAbs(p.addMod("epsAbs",mod,"ODE stepper absolute precision",epsAbsDefault)),
    sf(p.addMod("steppingFunction",mod,"Stepping function for EvolvedGSL",evolved::SF_RKCK)),
    nextDtTryCorretionFactor(p.addMod("nextDtTryCorretionFactor",mod,"Avoiding the tiny-next-timestep-to-try effect at the end of time intervals",100.))
{}


ParsStochastic::ParsStochastic(parameters::ParameterTable& p, const std::string& mod)
  : ParsEvolved(p,mod),
    seed(p.addTitle("StochasticTrajectory",mod).addMod("seed",mod,"Random number generator seed",1001ul)),
    noise(p.addMod("noise",mod,"Switching noise on/off",true)),
    nTraj(p.addMod("nTraj",mod,"Number of trajectories",size_t(100))) 
{}


} // trajectory