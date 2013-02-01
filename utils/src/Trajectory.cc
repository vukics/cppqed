#include "Trajectory.h"

#include "ParsTrajectory.h"
#include "impl/FormDouble.tcc"

#include <fstream>
#include <iostream>
#include <iomanip>


using namespace std;


namespace trajectory {


void run(Trajectory& traj, const Pars& p)
{ 
  if (!p.Dt) {cerr<<"Nonzero Dt required!"<<endl; return;}
  if (p.NDt) run(traj,p.NDt,p.Dt,p.displayInfo);
  else       run(traj,p.T  ,p.Dt,p.displayInfo);
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


Trajectory::Trajectory(const Pars& p)
  : ostream_(TrajectoryBaseHelper(p.ofn,p.precision)),
    precision_(p.precision)
{
} 


ostream& Trajectory::display(ostream& os, int precision) const
{
  const FormDouble fd(formdouble::positive(precision));
  return display_v( os<<fd(getTime())<<fd(getDtDid()) , precision)
           <<std::endl; // Note: endl flushes the buffer
}


ostream& Trajectory::displayKey(ostream& os) const
{
  return displayKey_v( os<<"# Trajectory\n#  1. time\n#  2. dtDid\n" );
}


namespace details {

void doRun(Trajectory& traj, double time, double deltaT)
{
  while (traj.getTime()<time) {
    traj.evolve(std::min(deltaT,time-traj.getTime()));
    traj.display();
  }
}

void doRun(Trajectory& traj, long nDt, double deltaT)
{
  for (long i=0; i<nDt; i++){
    traj.evolve(deltaT);
    traj.display();
  }
}

} // details

} // trajectory
