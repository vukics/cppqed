// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "Evolved.tcc"
#include "FormDouble.h"
#include "ParsTrajectory.h"
#include "SmartPtr.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <numeric> // for accumulate

namespace trajectory {


template<typename A>
void run(Adaptive<A>& traj, const ParsRun& p)
{
  if      (p.dc) run(traj,p.T,p.dc,p.sdf,p.ofn,p.initialFileName,p.precision,p.displayInfo,p.firstStateDisplay);
  else if (p.Dt) run(static_cast<Trajectory&>(traj),p);
  else std::cerr<<"Nonzero dc OR Dt required!"<<std::endl;
}


namespace details {

void writeNextArchive(std::ofstream*, const std::ostringstream&);
void readNextArchive(std::ifstream&, std::istringstream&);

} // details


template<typename T>
void writeViaSStream(const T& traj, std::ofstream* ofs)
{
  if (ofs && ofs->is_open()) {
    std::ostringstream oss(std::ios_base::binary);
    cpputils::oarchive stateArchive(oss);
    traj.writeState(stateArchive);
    details::writeNextArchive(ofs,oss);
  }
}

template<typename T>
void readViaSStream(T& traj, std::ifstream& ifs)
{
  std::istringstream iss(std::ios_base::binary);
  details::readNextArchive(ifs,iss);
  cpputils::iarchive stateArchive(iss);
  traj.readState(stateArchive);
}

namespace details {

  
namespace runTraits {

inline bool doContinue(const Trajectory& traj, double time, long      ) {return traj.getTime()<time;}
inline bool doContinue(const Trajectory&     , long length, long count) {return count<=length      ;}

inline void advance(Trajectory & traj, long       , double deltaT) {traj.evolve(deltaT);}
inline void advance(Trajectory & traj, double time, double deltaT) {traj.evolve(std::min(deltaT,time-traj.getTime()));}

template<typename A>
inline void advance(Adaptive<A>& traj, double time, int          ) {traj.step(time-traj.getTime());}

inline bool doDisplay(long      , double         ) {return true;}
inline bool doDisplay(long count, int displayFreq) {return !(count%displayFreq);}

inline const std::string writeTimestep(int   ) {return " timestep";}
inline const std::string writeTimestep(double) {return ""         ;}

inline double endTime(long    nDt, double dt, double currentTime) {return nDt*dt+currentTime;}
inline double endTime(double time, double   , double            ) {return time              ;}

} // runTraits


bool restoreState(Trajectory&, const std::string&, const std::string&, const std::string&);


template<typename T, typename L, typename D>
void run(T& traj, L length, D displayFreq, unsigned stateDisplayFreq, const std::string& trajectoryFileName, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay)
{
  using namespace std; using namespace boost; using namespace runTraits; using namespace cpputils;

  ////////////////////////////////////////////////
  // Determining i/o streams, eventual state input
  ////////////////////////////////////////////////

  static const string stateExtension(".state");
  const string stateFileName(trajectoryFileName+stateExtension);
  
  const bool
    outputToFile=(trajectoryFileName!=""),  
    continuing=restoreState(traj,trajectoryFileName,stateFileName,initialFileName);
  
  const double timeToReach=endTime(length,displayFreq,traj.getTime());
    
  if (timeToReach && timeToReach<=traj.getTime()) return;

  const boost::shared_ptr<ostream> outstream(!outputToFile ?
                                             nonOwningSharedPtr<ostream>(&cout) :
                                             static_pointer_cast<ostream>(make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app))); // regulates the deletion policy
  
  ostream& os=*outstream;

  if (os.fail()) throw TrajectoryFileOpeningException(trajectoryFileName);
 
  os<<setprecision(formdouble::actualPrecision(precision));
  
  ///////////////////////
  // Writing introduction
  ///////////////////////

  if (displayInfo) {
    if (!continuing)
      traj.displayParameters(os)
        <<endl<<"# Run Trajectory up to time "<<timeToReach
        <<" -- Display period: "<<displayFreq<<writeTimestep(displayFreq)<<endl<<endl;
    else
      os<<"# Continuing from time "<<traj.getTime()<<" up to time "<<timeToReach<<endl;
  }

  if (!timeToReach) {traj.display(os,precision); return;}

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  const boost::shared_ptr<ofstream> ofs = !outputToFile ? make_shared<ofstream>() : make_shared<ofstream>(stateFileName.c_str(),ios_base::app);
  bool stateSaved=false, evsDisplayed=false;

  try {

    for (long count=0, stateCount=0; doContinue(traj,length,count); ++count) {

      if (count) {
        advance(traj,length,displayFreq);
        stateSaved=evsDisplayed=false;
      }

      if (!count || doDisplay(count,displayFreq)) {

        if (
            stateDisplayFreq && 
            !(stateCount%stateDisplayFreq) && 
            (stateCount || (firstStateDisplay && !continuing))
           )
        {
          writeViaSStream(traj,ofs.get());
          stateSaved=true;
        }
        ++stateCount;

        if (count || !continuing) traj.display(os,precision);
        evsDisplayed=true;
      }
    }

  } catch (const StoppingCriterionReachedException& except) {os<<"\n# Stopping criterion has been reached"<<endl;}
  if (!evsDisplayed) traj.display(os,precision);

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  traj.logOnEnd(os);
  if (!stateSaved)   writeViaSStream(traj,ofs.get());
  
}

} // details


void run(Trajectory & traj, double time, double deltaT, unsigned sdf, const std::string& ofn, std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay)
{details::run(traj,time,deltaT,sdf,ofn,initialFileName,precision,displayInfo,firstStateDisplay);}

void run(Trajectory & traj, long   nDt , double deltaT, unsigned sdf, const std::string& ofn, std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay)
{details::run(traj,nDt ,deltaT,sdf,ofn,initialFileName,precision,displayInfo,firstStateDisplay);}

template<typename A>
void run(Adaptive<A>& traj, double time, int    dc    , unsigned sdf, const std::string& ofn, std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay)
{details::run(traj,time,dc    ,sdf,ofn,initialFileName,precision,displayInfo,firstStateDisplay);}


template<typename A>
AdaptiveIO<A>::AdaptiveIO(typename EvolvedIO::Ptr evolvedIO)
  : meta_(cpputils::TypeID<A>::value,
          SerializationMetadata::ARRAY_ONLY,
          cpputils::Rank<A>::value),
    evolvedIO_(evolvedIO)
{};

template<typename A>
cpputils::iarchive& AdaptiveIO<A>::readState(cpputils::iarchive& iar)
{
  bool dimension_check = meta_.trajectoryID != SerializationMetadata::ARRAY_ONLY;
  iar & meta_;
  if (meta_.rank!=cpputils::Rank<A>::value)
    throw RankMismatchException();
  std::vector<size_t> dims = cpputils::dimensions(evolvedIO_->getA());
  iar & *evolvedIO_;
  if (dimension_check && dims != cpputils::dimensions(evolvedIO_->getA()))
      throw DimensionsMismatchException();
  return iar;
}

template<typename A>
cpputils::oarchive& AdaptiveIO<A>::writeState(cpputils::oarchive& oar) const
{
  return oar & meta_ & *evolvedIO_;
}


template<typename A>
Adaptive<A>::Adaptive(A& y, typename Evolved::Derivs derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs, const evolved::Maker<A>& maker)
  : AdaptiveIO<A>(maker(y,derivs,dtInit,epsRel,epsAbs,scaleAbs)),
    evolved_(boost::dynamic_pointer_cast<Evolved>(AdaptiveIO<A>::getEvolvedIO())),
    dtInit_(dtInit)
{}


template<typename A>
Adaptive<A>::Adaptive(A& y, typename Evolved::Derivs derivs, double dtInit, const ParsEvolved& p, const A& scaleAbs, const evolved::Maker<A>& maker)
  : Adaptive(y,derivs,dtInit,p.epsRel,p.epsAbs,scaleAbs,maker) {}


template<typename A>
std::ostream& Adaptive<A>::displayParameters_v(std::ostream& os) const 
{
  return evolved_->displayParameters(os<<std::endl)<<"# Trajectory Parameters: epsRel="<<evolved_->getEpsRel()<<" epsAbs="<<evolved_->getEpsAbs()<<std::endl;
}

template<typename A>
cpputils::iarchive&  Adaptive<A>::readState_v(cpputils::iarchive& iar)
{
  AdaptiveIO<A>::readState(iar);
  if (meta_.trajectoryID != SerializationMetadata::ARRAY_ONLY) {
    if(meta_.trajectoryID != trajectoryID())
      throw TrajectoryMismatchException();
    readStateMore_v(iar);
  }
  if (getDtTry()==0)
    evolved_->setDtTry(dtInit_); // reset cached initial dtTry
  return iar;
}

template<typename A>
cpputils::oarchive&  Adaptive<A>::writeState_v(cpputils::oarchive& oar) const
{
  meta_.trajectoryID = trajectoryID(); // it is set here rather than @ construction as it is not good to call virtual functions @ construction
  AdaptiveIO<A>::writeState(oar);
  return writeStateMore_v(oar);
}

} // trajectory

#endif // CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED
