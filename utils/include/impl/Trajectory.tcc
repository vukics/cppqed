// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "impl/Evolved.tcc"
#include "ParsTrajectory.h"
#include "SmartPtr.h"

#include <boost/concept_check.hpp>

namespace trajectory {


template<typename A>
void run(Adaptive<A>& traj, const ParsRun& p)
{
  if      (p.dc) run(traj,p.T,p.dc,p.displayInfo);
  else if (p.Dt) run(static_cast<Trajectory&>(traj),p);
  else std::cerr<<"Nonzero dc OR Dt required!"<<std::endl;
}


namespace details {

  
namespace runTraits {

inline bool doContinue(const Trajectory& traj, double time, long      ) {return traj.getTime()<time;}
inline bool doContinue(const Trajectory&     , long length, long count) {return count<length       ;}

inline void advance(Trajectory & traj, long       , double deltaT) {traj.evolve(deltaT);}
inline void advance(Trajectory & traj, double time, double deltaT) {traj.evolve(std::min(deltaT,time-traj.getTime()));}

template<typename A>
inline void advance(Adaptive<A>& traj, double time, int          ) {traj.step(time-traj.getTime());}

inline bool doDisplay(long      , double         ) {return true;}
inline bool doDisplay(long count, int displayFreq) {return !(count%displayFreq);}

inline const std::string writeTimestep(int   ) {return " timestep";}
inline const std::string writeTimestep(double) {return ""         ;}

}

template<typename T, typename L, typename D>
void run(T& traj, L length, D displayFreq, const std::string& ofn, int precision, bool displayInfo)
{
  using namespace std; using namespace boost; using namespace runTraits;

  const bool continuing=false;
  
  struct {
    const shared_ptr<ostream> os_;
  } outstreams = { ofn=="" ? cpputils::nonOwningSharedPtr<ostream>(&cout) : shared_ptr<ostream>(new ofstream(ofn.c_str(),ios_base::app)) };
  
  ostream& os=*outstreams.os_;

  if (os.fail()) throw OutfileOpeningException(ofn);
 
  os<<setprecision(formdouble::actualPrecision(precision));
  
  if (displayInfo) {
    if (!continuing)
      traj.displayKey(traj.displayParameters(os)<<endl<<"# Key to data:"<<endl)<<endl<<"# Run Trajectory. Displaying in every "<<displayFreq<<writeTimestep(displayFreq)<<endl<<endl;
    else
      os<<"# Continuing..."<<endl;
  }

  try {
    if (!continuing) traj.display(os,precision);
    for (long count=0; doContinue(traj,length,count); ++count) {
      advance(traj,length,displayFreq);
      if (doDisplay(count,displayFreq)) traj.display(os,precision);
    }
  }
  catch (const StoppingCriterionReachedException& except) {
    os<<"# Stopping criterion has been reached"<<endl;
    throw except;
  }
}

} // details


void run(Trajectory & traj, double time, double deltaT, const std::string& ofn, int precision, bool displayInfo) {details::run(traj,time,deltaT,ofn,precision,displayInfo);}

void run(Trajectory & traj, long   nDt , double deltaT, const std::string& ofn, int precision, bool displayInfo) {details::run(traj,nDt ,deltaT,ofn,precision,displayInfo);}

template<typename A>
void run(Adaptive<A>& traj, double time, int    dc    , const std::string& ofn, int precision, bool displayInfo) {details::run(traj,time,dc    ,ofn,precision,displayInfo);}




template<typename A>
Adaptive<A>::Adaptive(A& y, typename Evolved::Derivs derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs, const evolved::Maker<A>& maker)
  : evolved_(maker(y,derivs,dtInit,epsRel,epsAbs,scaleAbs))
{}


template<typename A>
Adaptive<A>::Adaptive(A& y, typename Evolved::Derivs derivs, double dtInit, const ParsEvolved& p,         const A& scaleAbs, const evolved::Maker<A>& maker)
  : evolved_(maker(y,derivs,dtInit,p.epsRel,p.epsAbs,scaleAbs))
{}


template<typename A>
std::ostream& Adaptive<A>::displayParameters_v(std::ostream& os) const 
{
  return evolved_->displayParameters(os<<std::endl)<<"# Trajectory Parameters: epsRel="<<evolved_->getEpsRel()<<" epsAbs="<<evolved_->getEpsAbs()<<std::endl;
}



} // trajectory

#endif // UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED
