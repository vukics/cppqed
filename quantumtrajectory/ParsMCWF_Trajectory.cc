#include "ParsMCWF_Trajectory.h"

#include "impl/Pars.tcc"


namespace quantumtrajectory {


ParsMCWF_Trajectory::ParsMCWF_Trajectory(parameters::ParameterTable& p, const std::string& mod)
  : ParsStochasticTrajectory(p,mod),
    dpLimit(p.addTitle("MCWF_Trajectory",mod).addMod("dpLimit",mod,"MCWFS stepper total jump probability limit",0.1)),
    overshootTolerance(p.addMod("overshootTolerance",mod,"Jump probability overshoot tolerance factor",10.)),
    svdc(p.addMod("svdc",mod,"Number of displays between two state-vector Displays",0u)),
    firstSVDisplay(p.addMod("firstSVDisplay",mod,"Displays state vector at startup",true)),
    svdPrecision(p.addMod("svdPrecision",mod,"Precision of state-vector display (0 -- use the overall precision)",0)),
#ifdef USE_BOOST_SERIALIZATION
    binarySVFile(p.addMod("binarySVFile",mod,"Statevector file is binary (.svbin, full precision) or text (.sv, lower precision)",false)),
#endif // USE_BOOST_SERIALIZATION
    initFile (p.addMod<std::string>("initFile" ,mod,"file containing the initial state vector","")),
    basisDim(p.addMod("basisDim",mod,"number of basis vectors the stochastic wave function is compared against",size_t(0))),
    basisFile(p.addMod<std::string>("basisFile",mod,"file containing the basis vectors","")),
    logLevel(p.addMod("logLevel",mod,"logging level",0))
{}


} // quantumtrajectory
