#include "Trajectory.h"

#include "ParsTrajectory.h"
#include "impl/FormDouble.tcc"

#include <fstream>
#include <iostream>
#include <iomanip>


using namespace std;


namespace trajectory {


Trajectory::~Trajectory()
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


Trajectory::Trajectory(ostream& os, int precision)
  : ostream_(os<<setprecision(formdouble::actualPrecision(precision))),
    precision_(precision)
{
}


Trajectory::Trajectory(const string& ofn, int precision)
  : ostream_(TrajectoryBaseHelper(ofn,precision)),
    precision_(precision)
{
}


Trajectory::Trajectory(const ParsTrajectory& p)
  : ostream_(TrajectoryBaseHelper(p.ofn,p.precision)),
    precision_(p.precision)
{
} 


void Trajectory::display() const
{
  const FormDouble fd(formdouble::positive(precision_));
  getOstream()<<fd(getTime())<<fd(getDtDid());
  displayMore();
  getOstream().flush();
}


void Trajectory::displayKey() const
{
  getOstream()<<"# Trajectory\n#  1. time\n#  2. dtDid\n";
  displayMoreKey();
}



void details::doRun(Trajectory& traj, double time, double deltaT)
{
  while (traj.getTime()<time) {
    traj.evolve(std::min(deltaT,time-traj.getTime()));
    traj.display();
  }
}

void details::doRun(Trajectory& traj, long nDt, double deltaT)
{
  for (long i=0; i<nDt; i++){
    traj.evolve(deltaT);
    traj.display();
  }
}

} // trajectory
