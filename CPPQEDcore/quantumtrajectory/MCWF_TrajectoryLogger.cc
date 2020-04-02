// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MCWF_TrajectoryLogger.h"

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/numeric.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <iostream>
#include <cmath>


using namespace std;


quantumtrajectory::mcwf::Logger::Logger(int logLevel, bool isHamiltonian, size_t nLindblads)
  : logLevel_(logLevel), isHamiltonian_(isHamiltonian), nLindblads_(nLindblads),
    nSteps_(), nOvershot_(), nToleranceOvershot_(), nFailedSteps_(), nHamiltonianCalls_(),
    dpMaxOvershoot_(), dpToleranceMaxOvershoot_(), normMaxDeviation_(),
    traj_()
{}


ostream& quantumtrajectory::mcwf::Logger::onEnd(ostream& os) const
{
  if (logLevel_) {
    os<<"\n# Total number of steps: "<<nSteps_<<endl
      <<"\n# dpLimit overshot: "<<nOvershot_<<" times, maximal overshoot: "<<dpMaxOvershoot_
      <<"\n# dpTolerance overshot: "<<nToleranceOvershot_<<" times, maximal overshoot: "<<dpToleranceMaxOvershoot_<<endl
      <<"\n# Maximal deviation of norm from 1: "<<normMaxDeviation_<<endl;
    if (isHamiltonian_)
      os<<"\n# Failed ODE steps: "<<nFailedSteps_
        <<"\n# Number of Hamiltonian calls: "<<nHamiltonianCalls_<<endl;
    os<<"\n# MCWF Trajectory:\n";
      
    for (auto i : traj_) os<<"# "<<i.first<<"\t"<<i.second<<std::endl;
    // NEEDS_WORK the lambda expression of old doesn’t work with c++11
    // boost::for_each(traj_,os<<constant("# ")<<bind(&pair<double,size_t>::first,_1)<<constant("\t")<<bind(&pair<double,size_t>::second,_1)<<constant("\n"));
    os<<endl;
  }
  return os;
}


void quantumtrajectory::mcwf::Logger::step() {++nSteps_;}


void quantumtrajectory::mcwf::Logger::stepBack(double dp, double dtDid, double newDtTry, double t, bool logControl)
{
  ++nToleranceOvershot_;
  if (logControl) dpToleranceMaxOvershoot_=max(dpToleranceMaxOvershoot_,dp);
  if (logLevel_>2)
    cout<<"# dpTolerance overshot: "<<dp<<" stepping back to "<<t<<" timestep decreased: "<<dtDid<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::mcwf::Logger::overshot(double dp, double oldDtTry, double newDtTry, bool logControl)
{
  ++nOvershot_;
  if (logControl) dpMaxOvershoot_=max(dpMaxOvershoot_,dp);
  if (logLevel_>2)
    cout<<"# dpLimit overshot: "<<dp<<" timestep decreased: "<<oldDtTry<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::mcwf::Logger::processNorm(double norm)
{
  normMaxDeviation_=max(normMaxDeviation_,fabs(1-norm));
  // NEEDS_WORK this should be somehow weighed by the timestep
}


void quantumtrajectory::mcwf::Logger::jumpOccured(double t, size_t lindbladNo)
{
  traj_.push_back(make_pair(t,lindbladNo));
  if (logLevel_>1)
    cout<<"# Jump No. "<<lindbladNo<<" at time "<<t<<endl;
}


void quantumtrajectory::mcwf::Logger::logFailedSteps(size_t n)
{
  nFailedSteps_+=n;
  if (logLevel_>3)
    cout<<"# Number of failed steps in this timestep: "<<n<<endl;
}


void quantumtrajectory::mcwf::Logger::hamiltonianCalled()
{
  nHamiltonianCalls_++;
}



ostream& quantumtrajectory::ensemble::displayLog(ostream& os, const LoggerList& loggerList, size_t nBins, size_t nJumpsPerBin)
{
  using namespace boost;

#define AVERAGE_function(f) accumulate(loggerList | adaptors::transformed(bind(&Logger::f,_1)),0.)/loggerList.size()
#define MAX_function(f) max_element(loggerList, bind(std::less<double>(),bind(&Logger::f,_1),bind(&Logger::f,_2)) )->f

  os<<"\n# Average number of total steps: "<<AVERAGE_function(nSteps_)<<endl
    <<"\n# On average, dpLimit overshot: "<<AVERAGE_function(nOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpMaxOvershoot_)
    <<"\n# On average, dpTolerance overshot: "<<AVERAGE_function(nToleranceOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpToleranceMaxOvershoot_)<<endl
    <<"\n# Maximal deviation of norm from 1: "<<MAX_function(normMaxDeviation_)<<endl;
  if (loggerList.front().isHamiltonian_)
    os<<"\n# Average number of failed ODE steps: "<<AVERAGE_function(nFailedSteps_)
      <<"\n# Average number of Hamiltonian calls: "<<AVERAGE_function(nHamiltonianCalls_)<<endl;
      
#undef  MAX_function
#undef  AVERAGE_function
      
  using namespace accumulators;
  
  size_t nTotalJumps=accumulate(loggerList | adaptors::transformed( bind(&Logger::MCWF_Trajectory::size,bind(&Logger::traj_,_1)) ) ,0);

  // Heuristic: if nBins is not given (0), then the number of bins is determined such that the bins contain nJumpsPerBin samples on average
  if (!nBins && nTotalJumps<2*nJumpsPerBin) {
    os<<"\n# Too few jumps for a histogram\n";
  }
  else {
    size_t actualNumberOfBins=nBins ? nBins : nTotalJumps/nJumpsPerBin;
    accumulator_set<double, features<tag::density> > acc( tag::density::num_bins = actualNumberOfBins , tag::density::cache_size = nTotalJumps );

    // fill accumulator
    for (auto i : loggerList) for (auto j : i.traj_) acc(j.first);

    const iterator_range<std::vector<std::pair<double, double> >::iterator> histogram=density(acc);

    os<<"\n# Histogram of jumps. Number of bins="<<actualNumberOfBins+2<<endl;
    for (auto i : histogram)
      os<<"# "<<i.first<<"\t"<<i.second<<endl;
  }
 
  os<<"\n# Average number of total jumps: "<<nTotalJumps/double(loggerList.size())<<endl;
 
  return os;
}

