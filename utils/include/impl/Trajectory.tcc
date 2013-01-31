// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "impl/Evolved.tcc"
#include "ParsTrajectory.h"


namespace trajectory {


namespace details {


template<typename T, typename L, typename D>
void run(T& traj, L l, D d, void (*doRun)(T&,L,D), bool timestep, bool displayInfo)
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
    doRun(traj,l,d);
  }
  catch (const StoppingCriterionReachedException& except) {
    traj.getOstream()<<"# Stopping criterion has been reached"<<endl;
    throw except;
  }
}

void doRun(Trajectory&, long   nDt , double deltaT);
// Evolves the system on an Adaptive for nDt display intervals deltaT

void doRun(Trajectory&, double time, double deltaT);
// Evolves the system on an Adaptive up to time T and Displays in every deltaT

template<typename A> 
void doRun(Adaptive<A>& traj, double time, int dc)
{
  for (int count=1; traj.getTime()<time; count++) {
    traj.step(time-traj.getTime());
    if (!(count%dc)) traj.display();
  }
}

} // details

void runDt (Trajectory & traj, double time, double deltaT, bool displayInfo) {details::run(traj,time,deltaT,details::doRun,false,displayInfo);}

void runNDt(Trajectory & traj, long   nDt , double deltaT, bool displayInfo) {details::run(traj,nDt ,deltaT,details::doRun,false,displayInfo);}

template<typename A>
void run   (Adaptive<A>& traj, double time, int    dc    , bool displayInfo) {details::run(traj,time,dc    ,details::doRun,true ,displayInfo);}


template<typename A>
void evolve(Adaptive<A>& traj, const Pars& p)
{
  if      (p.dc) run  (traj,p.T,p.dc,p.displayInfo);
  else if (p.Dt) {
    if (p.NDt) runNDt(traj,p.NDt,p.Dt,p.displayInfo);
    else runDt(traj,p.T,p.Dt,p.displayInfo);
  }
  else std::cerr<<"Nonzero dc OR Dt required!"<<std::endl;
  //  std::cerr<<std::endl;
}


template<typename A>
Adaptive<A>::Adaptive(A& y, typename Evolved::Derivs derivs, double dtInit, double epsRel, double epsAbs,
                      const A& scaleAbs, const evolved::Maker<A>& maker)
  : evolved_(maker(y,derivs,dtInit,epsRel,epsAbs,scaleAbs))
{}


template<typename A>
Adaptive<A>::Adaptive(A& y, typename Evolved::Derivs derivs, double dtInit,
                      const A& scaleAbs, const Pars& p, const evolved::Maker<A>& maker)
  : evolved_(maker(y,derivs,dtInit,p.epsRel,p.epsAbs,scaleAbs)) {}


template<typename A>
void Adaptive<A>::displayParameters_v() const 
{
  evolved_->displayParameters(getOstream()<<std::endl)<<"# Trajectory Parameters: precision="<<getPrecision()<<" epsRel="<<evolved_->getEpsRel()<<" epsAbs="<<evolved_->getEpsAbs()<<std::endl;
}



} // trajectory

#endif // UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED
