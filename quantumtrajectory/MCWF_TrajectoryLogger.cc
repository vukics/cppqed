#include "MCWF_TrajectoryLogger.h"

#include "range_ex/algorithm.hpp"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <iostream>
#include <cmath>


using namespace std;


quantumtrajectory::MCWF_TrajectoryLogger::MCWF_TrajectoryLogger(unsigned logLevel, std::ostream& os)
  : logLevel_(logLevel), os_(os), nSteps_(), nOvershot_(), nToleranceOvershot_(), dpMaxOvershoot_(), dpToleranceMaxOvershoot_(), normMaxDeviation_(), traj_()
{}


quantumtrajectory::MCWF_TrajectoryLogger::~MCWF_TrajectoryLogger()
{
  using namespace boost::lambda;
  if (logLevel_) {
    os_
      <<  "# Total number of steps: "<<nSteps_
      <<"\n# dpLimit overshot: "<<nOvershot_<<" times, maximal overshoot: "<<dpMaxOvershoot_
      <<"\n# dpTolerance overshot: "<<nToleranceOvershot_<<" times, maximal overshoot: "<<dpToleranceMaxOvershoot_
      <<"\n# Maximal deviation of norm from 1: "<<normMaxDeviation_<<endl
      <<"# Trajectory:\n";
    boost::for_each(traj_,os_<<constant("# ")<<bind(&pair<double,size_t>::first,_1)<<constant("\t")<<bind(&pair<double,size_t>::second,_1)<<constant("\n"));
  }
}


void quantumtrajectory::MCWF_TrajectoryLogger::step() const {++nSteps_;}


void quantumtrajectory::MCWF_TrajectoryLogger::stepBack(double dp, double dtDid, double newDtTry, double t) const
{
  ++nToleranceOvershot_;
  dpToleranceMaxOvershoot_=max(dpToleranceMaxOvershoot_,dp);
  if (logLevel_>2)
    os_<<"# dpTolerance overshot: "<<dp<<" stepping back to "<<t<<" timestep decreased: "<<dtDid<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::MCWF_TrajectoryLogger::overshot(double dp, double oldDtTry, double newDtTry) const
{
  ++nOvershot_;
  dpMaxOvershoot_=max(dpMaxOvershoot_,dp);
  if (logLevel_>2)
    os_<<"# dpLimit overshot: "<<dp<<" timestep decreased: "<<oldDtTry<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::MCWF_TrajectoryLogger::processNorm(double norm) const
{
  normMaxDeviation_=max(normMaxDeviation_,fabs(1-norm));
  // NEEDS_WORK this should be somehow weighed by the timestep
}


void quantumtrajectory::MCWF_TrajectoryLogger::jumpOccured(double t, size_t jumpNo) const
{
  traj_.push_back(make_pair(t,jumpNo));
  if (logLevel_>1)
    os_<<"# Jump No. "<<jumpNo<<" at time "<<t<<endl;
}
