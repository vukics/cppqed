// -*- C++ -*-
#ifndef   _TRAJECTORY_IMPL_H
#define   _TRAJECTORY_IMPL_H

#include "ParsTrajectory.h"


namespace trajectory {


namespace details {


template<typename T, typename D>
void run(T& traj, double time, D d, void (*doRun)(T&,double,D), bool timestep, bool displayInfo)
{
  using namespace std;
  bool continuing=traj.getTime();
  
  if (displayInfo) {
    if (!continuing) {
      traj.displayParameters();
      traj.getOstream()<<endl<<"# Key to data:"<<endl;
      traj.displayKey();
      traj.getOstream()<<endl<<"# Run Trajectory. Displaying in every "<<d<<(timestep ? " timestep" : "")<<endl<<endl;
    }
    else {
      traj.getOstream()<<"# Continuing..."<<endl;
    }
  }

  try {
    if (!continuing) traj.display(); 
    doRun(traj,time,d);
  }
  catch (const StoppingCriterionReachedException& except) {
    traj.getOstream()<<"# Stopping criterion has been reached"<<endl;
    throw except;
  }
}

void doRun(TrajectoryBase&, double time, double deltaT);
// Evolves the system on a Trajectory up to time T and Displays in every deltaT

template<typename A> 
void doRun(Trajectory<A>& traj, double time, int dc)
{
  for (int count=1; traj.getTime()<time; count++) {
    traj.step(time-traj.getTime());
    if (!(count%dc)) traj.display();
  }
}

} // details

inline void runDt(TrajectoryBase& traj, double time, double deltaT, bool displayInfo) {details::run(traj,time,deltaT,details::doRun,false,displayInfo);}

template<typename A>
inline void run  (Trajectory<A> & traj, double time, int    dc    , bool displayInfo) {details::run(traj,time,dc    ,details::doRun,true ,displayInfo);}


template<typename A>
inline void evolve(Trajectory<A>& traj, const ParsTrajectory& p)
{
  if      (p.dc) run  (traj,p.T,p.dc,p.displayInfo);
  else if (p.Dt) runDt(traj,p.T,p.Dt,p.displayInfo);
  else std::cerr<<"Nonzero dc OR Dt required!"<<std::endl;
  std::cerr<<std::endl;
}


template<typename A>
Trajectory<A>::Trajectory(A& y, typename Evolved::Derivs derivs, double dtInit, double epsRel, double epsAbs,
			  const A& scaleAbs, const evolved::Maker<A>& maker)
  : evolved_(maker(y,derivs,dtInit,epsRel,epsAbs,scaleAbs))
{}


template<typename A>
Trajectory<A>::Trajectory(A& y, typename Evolved::Derivs derivs, double dtInit,
			  const A& scaleAbs, const ParsTrajectory& p, const evolved::Maker<A>& maker)
  : evolved_(maker(y,derivs,dtInit,p.epsRel,p.epsAbs,scaleAbs)) {}


template<typename A>
void Trajectory<A>::displayParameters() const 
{
  evolved_->displayParameters(getOstream()<<std::endl)<<"# Trajectory Parameters: epsRel="<<evolved_->getEpsRel()<<" epsAbs="<<evolved_->getEpsAbs()<<std::endl;
}



} // trajectory

#endif // _TRAJECTORY_IMPL_H
