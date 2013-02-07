// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "FormDouble.h"
#include "impl/Evolved.tcc"
#include "ParsTrajectory.h"
#include "SmartPtr.h"


namespace trajectory {


template<typename A>
void run(Adaptive<A>& traj, const ParsRun& p)
{
  if      (p.dc) run(traj,p.T,p.dc,p.ofn,p.precision,p.displayInfo);
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
inline bool doDisplay(long count, int displayFreq) {return !((count+1)%displayFreq);}

inline const std::string writeTimestep(int   ) {return " timestep";}
inline const std::string writeTimestep(double) {return ""         ;}

inline double endTime(long    nDt, double dt, double currentTime=0.) {return nDt*dt+currentTime;}
inline double endTime(double time, double   , double            =0.) {return time              ;}

}

template<typename T, typename L, typename D>
void run(T& traj, L length, D displayFreq, const std::string& trajectoryFileName, int precision, bool displayInfo)
{
  using namespace std; using namespace boost; using namespace runTraits;

  static const string stateExtension(".state");
  const string stateFileName(trajectoryFileName+stateExtension);
  
  bool continuing=false;
  
  if (trajectoryFileName!="") {
    ifstream trajectoryFile(trajectoryFileName.c_str());
    if (trajectoryFile.is_open() && trajectoryFile.peek()!=EOF) {
      ifstream stateFile(stateFileName.c_str());
      if (!stateFile.is_open()) throw StateFileOpeningException(stateFileName);
      cpputils::iarchive stateArchive(stateFile);
      traj.readState(stateArchive);
      if (endTime(length,displayFreq,traj.getTime())<=traj.getTime()) return;
      continuing=true;
    }
  }

  const shared_ptr<ostream> outstream(trajectoryFileName=="" ? cpputils::nonOwningSharedPtr<ostream>(&cout) : shared_ptr<ostream>(new ofstream(trajectoryFileName.c_str(),ios_base::app)));
  // regulates the deletion policy
  
  ostream& os=*outstream;

  if (os.fail()) throw TrajectoryFileOpeningException(trajectoryFileName);
 
  os<<setprecision(formdouble::actualPrecision(precision));
  
  if (displayInfo) {
    if (!continuing)
      traj.displayKey(traj.displayParameters(os)<<endl<<"# Key to data:"<<endl)
        <<endl<<"# Run Trajectory for time "<<endTime(length,displayFreq)
        <<" -- Display period: "<<displayFreq<<writeTimestep(displayFreq)<<endl<<endl;
    else
      os<<"# Continuing up to time "<<endTime(length,displayFreq,traj.getTime())<<endl;
  }

  try {
    if (!continuing) traj.display(os,precision);
    for (long count=0; doContinue(traj,length,count); ++count) {
      advance(traj,length,displayFreq);
      if (doDisplay(count,displayFreq)) traj.display(os,precision);
    }
#define ON_END \
    traj.logOnEnd(os); \
    if (trajectoryFileName!="") { \
      ofstream stateFile(stateFileName.c_str()); \
      cpputils::oarchive stateArchive(stateFile); \
      traj.writeState(stateArchive); \
    } 
    ON_END
  }
  catch (const StoppingCriterionReachedException& except) {
    os<<"# Stopping criterion has been reached"<<endl;
    ON_END
#undef ON_END
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
