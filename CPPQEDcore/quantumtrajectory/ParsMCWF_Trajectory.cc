// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsMCWF_Trajectory.h"

#include "Pars.h"


quantumtrajectory::mcwf::Pars::Pars(parameters::Table& p, const std::string& mod)
  : ParsStochastic(p,mod),
    dpLimit(p.addTitle("MCWF_Trajectory",mod).addMod("dpLimit",mod,"MCWFS stepper total jump probability limit",0.01)),
    overshootTolerance(p.addMod("overshootTolerance",mod,"Jump probability overshoot tolerance factor",10.)),
    nBins(p.addMod("nBins",mod,"number of bins used for the histogram of jumps created by EnsembleMCWF",size_t(0))),
    nJumpsPerBin(p.addMod("nJumpsPerBin",mod,"average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination",size_t(50)))
{}

