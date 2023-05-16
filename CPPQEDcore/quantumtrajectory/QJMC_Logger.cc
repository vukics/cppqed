// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "QJMC_Logger.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <cmath>
#include <iostream>
#include <numeric>

using namespace std;

using cppqedutils::LogTree;

quantumtrajectory::qjmc::Logger::Logger(size_t nLindblads) : nLindblads_(nLindblads) {}


LogTree quantumtrajectory::qjmc::Logger::outro() const
{
  return {
    {"Total number of MCWF steps",nQJMC_steps_},
    {"dpLimit overshot",{
      {"nOvershot",nOvershot_},
      {"maximal overshoot",dpMaxOvershoot_}}},
    {"normMaxDeviation",normMaxDeviation_},
    {"Jump trajectory",traj_}};
}


void quantumtrajectory::qjmc::Logger::step() {++nQJMC_steps_;}


LogTree quantumtrajectory::qjmc::Logger::overshot(double dp, double oldDtTry, double newDtTry)
{
  ++nOvershot_; dpMaxOvershoot_=max(dpMaxOvershoot_,dp);
  return {{"dpLimit overshot",dp},{"timestep decreased",{{"from",oldDtTry},{"to",newDtTry}}}};
}


void quantumtrajectory::qjmc::Logger::processNorm(double norm)
{
  normMaxDeviation_=max(normMaxDeviation_,fabs(1-norm));
  // TODO: this should be somehow weighed by the timestep
}


LogTree quantumtrajectory::qjmc::Logger::jumpOccured(double t, size_t lindbladNo)
{
  traj_.push_back({t,lindbladNo});
  return {{"no.",lindbladNo},{"at time",t}};
}


/*
ostream& quantumtrajectory::qjmc::EnsembleLogger::stream(ostream& os, const LoggerList& loggerList, size_t nBins, size_t nJumpsPerBin)
{
  using namespace boost;

#define AVERAGE_function(f) accumulate(loggerList.begin(),loggerList.end(), 0., [] (double init, const Logger& l) {return init + l.f;})/loggerList.size()
#define MAX_function(f) ranges::max_element(loggerList, [] (const Logger& a, const Logger& b) {return a.f<b.f;})->f

  os<<"\nAverage number of total steps: "<<AVERAGE_function(nMCWF_steps_)<<endl
    <<"\nOn average, dpLimit overshot: "<<AVERAGE_function(nOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpMaxOvershoot_)
    <<"\nOn average, dpTolerance overshot: "<<AVERAGE_function(nToleranceOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpToleranceMaxOvershoot_)<<endl
    <<"\nMaximal deviation of norm from 1: "<<MAX_function(normMaxDeviation_)<<endl;
// NEEDS_WORK: this was before factoring out logging to Evolved/Adaptive, so it could be restored on some more fundamental level:  if (loggerList.front().isHamiltonian_)
//    os<<"\nAverage number of failed ODE steps: "<<AVERAGE_function(nFailedSteps_)<<"\nAverage number of Hamiltonian calls:"<<AVERAGE_function(nHamiltonianCalls_)<<endl;
      
#undef  MAX_function
#undef  AVERAGE_function
      
  using namespace accumulators;
  
  size_t nTotalJumps=accumulate(loggerList.begin(), loggerList.end(), 0, [] (double init, const Logger& l) {return init + l.traj_.size();});

  // Heuristic: if nBins is not given (0), then the number of bins is determined such that the bins contain nJumpsPerBin samples on average
  if (!nBins && nTotalJumps<2*nJumpsPerBin) {
    os<<"\nToo few jumps for a histogram\n";
  }
  else {
    size_t actualNumberOfBins=nBins ? nBins : nTotalJumps/nJumpsPerBin;
    accumulator_set<double, features<tag::density> > acc( tag::density::num_bins = actualNumberOfBins , tag::density::cache_size = nTotalJumps );

    // fill accumulator
    for (auto i : loggerList) for (auto j : i.traj_) acc(j.first);

    const iterator_range<std::vector<std::pair<double, double> >::iterator> histogram=density(acc);

    os<<"\nHistogram of jumps. Number of bins="<<actualNumberOfBins+2<<endl;
    for (auto i : histogram)
      os<<i.first<<"\t"<<i.second<<endl;
  }
 
  os<<"\nAverage number of total jumps: "<<nTotalJumps/double(loggerList.size())<<endl;
 
  return os;
}

*/
