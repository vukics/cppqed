#include "Trajectory.h"

#include "ParsTrajectory.h"
#include "FormDouble.h"

#include <fstream>
#include <iostream>
#include <iomanip>


using namespace std;


namespace trajectory {


TrajectoryBase::~TrajectoryBase()
{
  if (dynamic_cast<ofstream*>(&ostream_)) delete &ostream_;
  // NEEDS_WORK
}


ostream& TrajectoryBaseHelper(const string& ofn, int precision)
{
  ofstream*const res = ( ofn=="" ? static_cast<ofstream*const>(&cout) : new ofstream(ofn.c_str(),ios_base::app) );

  if (res->fail()) throw OutfileOpeningException(ofn);

  return (*res)<<setprecision(formdouble::actualPrecision(precision));
}


TrajectoryBase::TrajectoryBase(ostream& os, int precision)
  : ostream_(os<<setprecision(formdouble::actualPrecision(precision))),
    precision_(precision)
{
}


TrajectoryBase::TrajectoryBase(const string& ofn, int precision)
  : ostream_(TrajectoryBaseHelper(ofn,precision)),
    precision_(precision)
{
}


TrajectoryBase::TrajectoryBase(const ParsTrajectory& p)
  : ostream_(TrajectoryBaseHelper(p.ofn,p.precision)),
    precision_(p.precision)
{
} 


void TrajectoryBase::display() const
{
  const FormDouble fd(formdouble::positive(precision_));
  getOstream()<<fd(getTime())<<fd(getDtDid());
  displayMore();
  getOstream().flush();
}


void TrajectoryBase::displayKey() const
{
  getOstream()<<"# Trajectory\n#  1. time\n#  2. dtDid\n";
  displayMoreKey();
}



void details::doRun(TrajectoryBase& traj, double time, double deltaT)
{
  while (traj.getTime()<time) {
    traj.evolve(std::min(deltaT,time-traj.getTime()));
    traj.display();
  }
}

void details::doRun(TrajectoryBase& traj, long nDt, double deltaT)
{
  for (long i=0; i<nDt; i++){
    traj.evolve(deltaT);
    traj.display();
  }
}

} // trajectory
