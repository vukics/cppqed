#include<fstream>
#include<iostream>

#include "ParsTrajectory.h"
#include "Trajectory.h"
#include "FormDouble.h"


using namespace std;


namespace trajectory {


TrajectoryBase::~TrajectoryBase()
{
  if (dynamic_cast<ofstream*>(&ostream_)) delete &ostream_;
  // NEEDS_WORK
}


ostream& TrajectoryBaseHelper(const string& ofn)
{
  if (ofn=="") return cout;
  ofstream*const res=new ofstream(ofn.c_str(),ios_base::app);
  if (res->fail()) throw OutfileOpeningException(ofn);
  return *res;
}


TrajectoryBase::TrajectoryBase(ostream& os, int precision)
  : ostream_(os),
    precision_(precision)
{
}


TrajectoryBase::TrajectoryBase(const string& ofn, int precision)
  : ostream_(TrajectoryBaseHelper(ofn)),
    precision_(precision)
{
}


TrajectoryBase::TrajectoryBase(const ParsTrajectory& p)
  : ostream_(TrajectoryBaseHelper(p.ofn)),
    precision_(p.precision)
{
} 


void TrajectoryBase::display() const
{
  const FormDouble fd(formdouble::positive(precision_));
  getOstream()<<fd(getTime())<<fd(getDtDid());
  displayMore(precision_);
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



} // trajectory
