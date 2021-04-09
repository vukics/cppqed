// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MCWF_TrajectoryLogger.h"

#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/numeric.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <iostream>
#include <cmath>


using namespace std;


quantumtrajectory::mcwf::Logger::Logger(int logLevel, size_t nLindblads)
  : logLevel_(logLevel), nLindblads_(nLindblads),
    nMCWF_steps_(), nOvershot_(), nToleranceOvershot_(),
    dpMaxOvershoot_(), dpToleranceMaxOvershoot_(), normMaxDeviation_(),
    traj_()
{}


ostream& quantumtrajectory::mcwf::Logger::onEnd(ostream& os) const
{
  if (logLevel_) {
    os<<"\nTotal number of MCWF steps: "<<nMCWF_steps_<<endl
      <<"\ndpLimit overshot: "<<nOvershot_<<" times, maximal overshoot: "<<dpMaxOvershoot_
      <<"\ndpTolerance overshot: "<<nToleranceOvershot_<<" times, maximal overshoot: "<<dpToleranceMaxOvershoot_<<endl
      <<"\nMaximal deviation of norm from 1: "<<normMaxDeviation_<<endl
      <<"\nMCWF Trajectory:\n";
      
    for (auto i : traj_) os<<i.first<<"\t"<<i.second<<std::endl;
    os<<endl;
  }
  return os;
}


void quantumtrajectory::mcwf::Logger::step() {++nMCWF_steps_;}


void quantumtrajectory::mcwf::Logger::stepBack(ostream& os, double dp, double dtDid, double newDtTry, double t, bool logControl)
{
  ++nToleranceOvershot_;
  if (logControl) dpToleranceMaxOvershoot_=max(dpToleranceMaxOvershoot_,dp);
  if (logLevel_>2)
    os<<"dpTolerance overshot: "<<dp<<" stepping back to "<<t<<" timestep decreased: "<<dtDid<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::mcwf::Logger::overshot(ostream& os, double dp, double oldDtTry, double newDtTry, bool logControl)
{
  ++nOvershot_;
  if (logControl) dpMaxOvershoot_=max(dpMaxOvershoot_,dp);
  if (logLevel_>2)
    os<<"dpLimit overshot: "<<dp<<" timestep decreased: "<<oldDtTry<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::mcwf::Logger::processNorm(double norm)
{
  normMaxDeviation_=max(normMaxDeviation_,fabs(1-norm));
  // NEEDS_WORK this should be somehow weighed by the timestep
}


void quantumtrajectory::mcwf::Logger::jumpOccured(ostream& os, double t, size_t lindbladNo)
{
  traj_.push_back({t,lindbladNo});
  if (logLevel_>1)
    os<<"Jump No. "<<lindbladNo<<" at time "<<t<<endl;
}



ostream& quantumtrajectory::mcwf::EnsembleLogger::stream(ostream& os, const LoggerList& loggerList, size_t nBins, size_t nJumpsPerBin)
{
  using namespace boost;

#define AVERAGE_function(f) accumulate(loggerList, 0., [] (double init, const Logger& l) {return init + l.f;})/loggerList.size()
#define MAX_function(f) max_element(loggerList, [] (const Logger& a, const Logger& b) {return a.f<b.f;})->f

  os<<"\nAverage number of total steps: "<<AVERAGE_function(nMCWF_steps_)<<endl
    <<"\nOn average, dpLimit overshot: "<<AVERAGE_function(nOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpMaxOvershoot_)
    <<"\nOn average, dpTolerance overshot: "<<AVERAGE_function(nToleranceOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpToleranceMaxOvershoot_)<<endl
    <<"\nMaximal deviation of norm from 1: "<<MAX_function(normMaxDeviation_)<<endl;
/* NEEDS_WORK: this was before factoring out logging to Evolved/Adaptive, so it could be restored on some more fundamental level:  if (loggerList.front().isHamiltonian_) os<<"\nAverage number of failed ODE steps: "<<AVERAGE_function(nFailedSteps_)<<"\nAverage number of Hamiltonian calls:"<<AVERAGE_function(nHamiltonianCalls_)<<endl;*/
      
#undef  MAX_function
#undef  AVERAGE_function
      
  using namespace accumulators;
  
  size_t nTotalJumps=accumulate(loggerList, 0, [] (double init, const Logger& l) {return init + l.traj_.size();});

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

