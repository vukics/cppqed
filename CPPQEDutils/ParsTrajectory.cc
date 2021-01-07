// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsStochasticTrajectory.h"

#include "FormDouble.h"
#include "Pars.h"

namespace trajectory {


const double ParsEvolved::epsRelDefault=1e-6;
const double ParsEvolved::epsAbsDefault=1e-12;


ParsRun::ParsRun(parameters::Table& p, const std::string& mod)
  : T(p.addTitle("Trajectory",mod).add("T",mod,"Simulated time",1.)),
    dc(p.add("dc",mod,"Number of steps between two streamings",10)),
    Dt(p.add("Dt",mod,"Timestep between two streamings",.1)),
    NDt(p.add("NDt",mod,"Number of steps in Dt mode",0L)),
    ofn(p.add<std::string>("o",mod,"Output file name for Trajectory, when empty, cout","")),
    initialFileName(p.add<std::string>("initialFileName",mod,"Trajectory initial file name","")),
    precision(p.add("precision",mod,"General precision of output",formdouble::Zero(FormDouble::defaultPrecision))),
    streamInfo(p.add("streamInfo",mod,"Whether to stream header for trajectories",true)),
    firstStateStream(p.add("firstStateStream",mod,"Streams trajectory state at startup",true)),
    sdf(p.add("sdf",mod,"State output frequency",0u)),
    autoStopEpsilon(p.add("autoStopEpsilon",mod,"Relative precision for autostopping",ParsEvolved::epsRelDefault)),
    autoStopRepetition(p.add("autoStopRepetition",mod,"Number of streamed lines repeated within relative precision before autostopping",0u)),
    parsedCommandLine_(p.getParsedCommandLine())
{}


ParsEvolved::ParsEvolved(parameters::Table& p, const std::string& mod)
  : epsRel(p.addTitle("Evolved",mod).add("eps"   ,mod,"ODE stepper relative precision",epsRelDefault)),
    epsAbs(p.add("epsAbs",mod,"ODE stepper absolute precision",epsAbsDefault)),
    sf(p.add("steppingFunction",mod,"Stepping function for EvolvedGSL",evolved::SF_RKCK)),
    nextDtTryCorrectionFactor(p.add("nextDtTryCorrectionFactor",mod,"Avoiding the tiny-next-timestep-to-try effect at the end of time intervals",100.)),
    logLevel(p.add("logLevel",mod,"logging level",0))
{}


ParsStochastic::ParsStochastic(parameters::Table& p, const std::string& mod)
  : ParsEvolved(p,mod),
    seed(p.addTitle("StochasticTrajectory",mod).add("seed",mod,"Random number generator seed",1001ul)),
    prngStream(p.add("prngStream",mod,"Random number generator independent stream ordinal",1ul)),
    noise(p.add("noise",mod,"Switching noise on/off",true)),
    nTraj(p.add("nTraj",mod,"Number of trajectories",size_t(100))) 
{}


} // trajectory
