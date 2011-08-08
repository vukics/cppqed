#include "MCWF_TrajectoryLogger.h"

#include <iostream>
#include <cmath>


using namespace std;


quantumtrajectory::MCWF_TrajectoryLogger::MCWF_TrajectoryLogger(unsigned logLevel, std::ostream& os)
  : logLevel_(logLevel), os_(os), nSteps_(), nOvershot_(), dpMaxOvershoot_(), normMaxDeviation_()
{}


quantumtrajectory::MCWF_TrajectoryLogger::~MCWF_TrajectoryLogger()
{
  if (logLevel_) 
    os_<<"# Total number of steps: "<<nSteps_<<"\n# dpLimit overshot: "<<nOvershot_<<" times, maximal overshoot: "<<dpMaxOvershoot_<<"\n# Maximal deviation of norm from 1: "<<normMaxDeviation_<<endl;
}


void quantumtrajectory::MCWF_TrajectoryLogger::step() const {++nSteps_;}


void quantumtrajectory::MCWF_TrajectoryLogger::overshot(double dp, double oldDtTry, double newDtTry) const
{
  ++nOvershot_;
  dpMaxOvershoot_=max(dpMaxOvershoot_,dp);
  if (logLevel_>1)
    os_<<"# dpLimit overshot: "<<dp<<" timestep decreased: "<<oldDtTry<<" => "<<newDtTry<<endl;
}


void quantumtrajectory::MCWF_TrajectoryLogger::processNorm(double norm) const
{
  normMaxDeviation_=max(normMaxDeviation_,fabs(1-norm));
  // NEEDS_WORK this should be somehow weighed by the timestep
}
