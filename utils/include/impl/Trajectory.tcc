// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "FormDouble.h"
#include "impl/Evolved.tcc"
#include "ParsTrajectory.h"
#include "SmartPtr.h"

#include <boost/make_shared.hpp>

#include <iomanip>
#include <fstream>

namespace trajectory {


template<typename A>
void run(Adaptive<A>& traj, const ParsRun& p)
{
  if      (p.dc) run(traj,p.T,p.dc,p.sdf,p.ofn,p.precision,p.displayInfo);
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
void run(T& traj, L length, D displayFreq, unsigned stateDisplayFreq, const std::string& trajectoryFileName, int precision, bool displayInfo)
{
  using namespace std; using namespace boost; using namespace runTraits; using namespace cpputils;

  ////////////////////////////////////////////////
  // Determining i/o streams, eventual state input
  ////////////////////////////////////////////////

  static const string stateExtension(".state");
  const string stateFileName(trajectoryFileName+stateExtension);
  
  const bool outputToFile=(trajectoryFileName!="");
  bool continuing=false;
  
  if (outputToFile) {
    ifstream trajectoryFile(trajectoryFileName.c_str());
    if (trajectoryFile.is_open() && (trajectoryFile.peek(), !trajectoryFile.eof()) ) {
      ifstream stateFile(stateFileName.c_str(), ios_base::binary); stateFile.exceptions(ifstream::eofbit);
      if (!stateFile.is_open()) throw StateFileOpeningException(stateFileName);
      iarchive stateArchive(stateFile);
      try {while (true) traj.readState(stateArchive);} catch (ifstream::failure) {} // boost::archive::archive_exception
      if (endTime(length,displayFreq,traj.getTime())<=traj.getTime()) return;
      continuing=true;
    }
  }

  const shared_ptr<ostream> outstream(!outputToFile ?
                                      nonOwningSharedPtr<ostream>(&cout) :
                                      static_pointer_cast<ostream>(make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app)));
  // regulates the deletion policy
  
  ostream& os=*outstream;

  if (os.fail()) throw TrajectoryFileOpeningException(trajectoryFileName);
 
  os<<setprecision(formdouble::actualPrecision(precision));
  
  ///////////////////////
  // Writing introduction
  ///////////////////////

  if (displayInfo) {
    if (!continuing)
      traj.displayKey(traj.displayParameters(os)<<endl<<"# Key to data:"<<endl)
        <<endl<<"# Run Trajectory up to time "<<endTime(length,displayFreq)
        <<" -- Display period: "<<displayFreq<<writeTimestep(displayFreq)<<endl<<endl;
    else
      os<<"# Continuing from time "<<traj.getTime()<<" up to time "<<endTime(length,displayFreq,traj.getTime())<<endl;
  }
  if (!continuing) traj.display(os,precision);

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  try {
    struct {
      const shared_ptr<ofstream> ofs_;
      const shared_ptr<oarchive> oar_;
    } oarchiveWrapper = {
                          (!stateDisplayFreq || !outputToFile) ? make_shared<ofstream>() : make_shared<ofstream>(stateFileName.c_str(), ios::binary | ios_base::app) ,
                          (!oarchiveWrapper.ofs_ || !oarchiveWrapper.ofs_->is_open()) ? shared_ptr<oarchive>() : make_shared<oarchive>(boost::ref(*oarchiveWrapper.ofs_))
                        };
      
    for (long count=0, stateCount=1; doContinue(traj,length,count); ++count) {
      advance(traj,length,displayFreq);
      if (doDisplay(count,displayFreq)) {traj.display(os,precision); ++stateCount;}
      if (oarchiveWrapper.oar_ && !(stateCount%stateDisplayFreq)) traj.writeState(*oarchiveWrapper.oar_);
    }
  } catch (const StoppingCriterionReachedException& except) {os<<"\n# Stopping criterion has been reached"<<endl;}

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  traj.logOnEnd(os);
  if (outputToFile) {
    ofstream stateFile(stateFileName.c_str(),ios::binary | ios_base::app);
    oarchive stateArchive(stateFile);
    traj.writeState(stateArchive);
  }
  
}

} // details


void run(Trajectory & traj, double time, double deltaT, unsigned sdf, const std::string& ofn, int precision, bool displayInfo) {details::run(traj,time,deltaT,sdf,ofn,precision,displayInfo);}

void run(Trajectory & traj, long   nDt , double deltaT, unsigned sdf, const std::string& ofn, int precision, bool displayInfo) {details::run(traj,nDt ,deltaT,sdf,ofn,precision,displayInfo);}

template<typename A>
void run(Adaptive<A>& traj, double time, int    dc    , unsigned sdf, const std::string& ofn, int precision, bool displayInfo) {details::run(traj,time,dc    ,sdf,ofn,precision,displayInfo);}




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
