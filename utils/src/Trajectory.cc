#include "Trajectory.h"

#include "ParsTrajectory.h"
#include "impl/FormDouble.tcc"

#include <fstream>
#include <iostream>
#include <iomanip>


using namespace std;


namespace trajectory {


void run(Trajectory& traj, const ParsRun& p)
{ 
  if (!p.Dt) {cerr<<"Nonzero Dt required!"<<endl; return;}
  if (p.NDt) run(traj,p.NDt,p.Dt,p.displayInfo);
  else       run(traj,p.T  ,p.Dt,p.displayInfo);
}


ostream& Trajectory::display(ostream& os, int precision) const
{
  const FormDouble fd(formdouble::positive(precision));
  return display_v( os<<fd(getTime())<<fd(getDtDid()) , precision)
           <<std::endl; // Note: endl flushes the buffer
}


ostream& Trajectory::displayKey(ostream& os) const
{
  size_t i=3;
  return displayKey_v( os<<"# Trajectory\n#  1. time\n#  2. dtDid\n" , i);
}


} // trajectory
