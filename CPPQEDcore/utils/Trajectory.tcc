// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "CommentingStream.h"
#include "Evolved.tcc"
#include "FormDouble.h"
#include "IO_Manip.h"
#include "ParsTrajectory.h"
#include "Version.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <numeric> // for accumulate


template<typename A, typename BASE>
void trajectory::run(Adaptive<A,BASE>& traj, const ParsRun& p)
{
  if      (p.dc) run(traj,p.T,p.dc,p.sdf,p.ofn,p.initialFileName,p.precision,p.displayInfo,p.firstStateDisplay,p.autoStopEpsilon,p.autoStopRepetition,p.getParsedCommandLine());
  else if (p.Dt) run(static_cast<Trajectory&>(traj),p);
  else std::cerr<<"Nonzero dc OR Dt required!"<<std::endl;
}


namespace trajectory { namespace details {

void writeNextArchive(std::ostream*, const std::ostringstream&);
void readNextArchive(std::istream*, std::istringstream&);

} // details


template<typename T>
void writeViaSStream(const T& traj, std::ostream* ofs)
{
  if (ofs) {
    std::ostringstream oss(std::ios_base::binary);
    cpputils::oarchive stateArchive(oss);
    traj.writeState(stateArchive);
    details::writeNextArchive(ofs,oss);
  }
}

template<typename T>
void readViaSStream(T& traj, std::istream* ifs)
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

template<typename A, typename BASE>
inline void advance(Adaptive<A,BASE>& traj, double time, int) {traj.step(time-traj.getTime());}

inline bool doDisplay(long      , double         ) {return true;}
inline bool doDisplay(long count, int displayFreq) {return !(count%displayFreq);}

inline const std::string writeTimestep(int   ) {return " timestep";}
inline const std::string writeTimestep(double) {return ""         ;}

inline double endTime(long    nDt, double dt, double currentTime) {return nDt*dt+currentTime;}
inline double endTime(double time, double   , double            ) {return time              ;}

} // runTraits


class DisplayAndAutostopHandler
{
public:
  typedef std::shared_ptr<const DisplayAndAutostopHandler> Ptr;

  DisplayAndAutostopHandler(const Trajectory& traj) : traj_(traj) {}

  virtual std::ostream& display(std::ostream& os, int precision) const; // const {return traj_.display(os);}

  virtual ~DisplayAndAutostopHandler() {}
  
protected:
  const Trajectory& traj_;

};

/** \todo The present design doesn’t provide real modularity, for this we would need a maker class for DisplayAndAutostopHandler. */
const DisplayAndAutostopHandler::Ptr makeDisplayAndAutostopHandler(const Trajectory&, double autoStopEpsilon, unsigned autoStopRepetition);

bool restoreState(Trajectory&, const std::string&, const std::string&, const std::string&);


template<typename T>
inline void nullDelete (T*) {}


template<typename T, typename L, typename D>
void run(T& traj, L length, D displayFreq, unsigned stateDisplayFreq, const std::string& trajectoryFileName, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay,
         double autoStopEpsilon, unsigned autoStopRepetition, const std::string& parsedCommandLine)
{
  using namespace std; using namespace runTraits; using namespace cpputils;

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

  const std::shared_ptr<ostream> outstream(!outputToFile ?
                                           std::shared_ptr<ostream>(&cout,details::nullDelete<ostream>) : // since cout is a system-wide object, this should be safe
                                           static_pointer_cast<ostream>(std::make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app))); // regulates the deletion policy
  
  if (outstream->fail()) throw TrajectoryFileOpeningException(trajectoryFileName);
  traj.setLogStreamDuringRun(outstream);
  
  ostream& os=*outstream;
  IO_Manipulator::_(os);
  os<<setprecision(formdouble::actualPrecision(precision));
  
  ///////////////////////
  // Writing introduction
  ///////////////////////

  {
  
  CommentingStream commentingStream(os);
  commentingStream<<setprecision(formdouble::actualPrecision(precision));
  
  // auto & commentingStream=os;
    
  if (displayInfo) {
    if (!continuing) {
      if (parsedCommandLine!="") commentingStream<<parsedCommandLine<<endl<<endl;
      traj.displayParameters(commentingStream<<versionHelper())
        <<endl<<"Run Trajectory up to time "<<timeToReach
        <<" -- Display period: "<<displayFreq<<writeTimestep(displayFreq)<<endl<<endl;
    }
    else
      commentingStream<<"Continuing from time "<<traj.getTime()<<" up to time "<<timeToReach<<endl;
  }

  }
  
  if (!timeToReach) {traj.display(os,precision); return;}

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  const std::shared_ptr<ostream> ofs = !outputToFile ? std::make_shared<ofstream>() : openStateFileWriting(stateFileName);

  bool
    stateSaved=false,   // signifies whether the state has already been saved for the actual time instant of the trajectory
    evsDisplayed=false; // signifies whether the expectation values have already been displayed ”

  auto startInstant=chrono::steady_clock::now();
  auto endInstant=chrono::steady_clock::now();
  
  try {

    auto displayAndAutostopHandler(makeDisplayAndAutostopHandler(traj,autoStopEpsilon,autoStopRepetition));

    startInstant=chrono::steady_clock::now();
    
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

        if (count || !continuing) {
          evsDisplayed=true;
          displayAndAutostopHandler->display(os,precision);
        }
      }
    }
    
    endInstant=chrono::steady_clock::now();

  } catch (const StoppingCriterionReachedException& except) {
    if (endInstant<startInstant) endInstant=chrono::steady_clock::now();
    CommentingStream(os)<<"Stopping criterion has been reached"<<endl;    
  }
  if (!evsDisplayed) traj.display(os,precision);

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  {
    CommentingStream temp(os);
    temp<<setprecision(formdouble::actualPrecision(precision))
        <<"Duration of main loop (us): "<<std::chrono::duration_cast<std::chrono::microseconds>(endInstant-startInstant).count()<<std::endl;
    traj.logOnEnd(temp);}
  if (!stateSaved) writeViaSStream(traj,ofs.get());
  
}

} } // trajectory::details


template<typename A, typename BASE>
void trajectory::run(Adaptive<A,BASE>& traj, double time, int dc, unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision,
                     bool displayInfo, bool firstStateDisplay,
                     double autoStopEpsilon, unsigned autoStopRepetition, const std::string& parsedCommandLine)
{details::run(traj,time,dc,sdf,ofn,initialFileName,precision,displayInfo,firstStateDisplay,autoStopEpsilon,autoStopRepetition,parsedCommandLine);}


template<typename A>
trajectory::AdaptiveIO<A>::AdaptiveIO(typename EvolvedIO::Ptr evolvedIO)
  : meta_(cpputils::TypeID<A>::value,
          SerializationMetadata::ARRAY_ONLY,
          cpputils::Rank<A>::value),
    evolvedIO_(evolvedIO)
{};

template<typename A>
cpputils::iarchive& trajectory::AdaptiveIO<A>::readState(cpputils::iarchive& iar)
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
cpputils::oarchive& trajectory::AdaptiveIO<A>::writeState(cpputils::oarchive& oar) const
{
  return oar & meta_ & *evolvedIO_;
}


template<typename A, typename BASE> template<typename... BaseInitializationPack>
trajectory::Adaptive<A,BASE>::Adaptive(A& y, Derivs derivs, double dtInit, int logLevel, double epsRel, double epsAbs, const A& scaleAbs, const evolved::Maker<A>& maker,
                                       BaseInitializationPack&&... bip)
  : AdaptiveIO<A>(maker(y,derivs,dtInit,epsRel,epsAbs,scaleAbs)),
    BASE(std::forward<BaseInitializationPack>(bip)...),
    evolved_(std::dynamic_pointer_cast<Evolved>(AdaptiveIO<A>::getEvolvedIO())),
    dtInit_(dtInit), logLevel_(logLevel)
{}


template<typename A, typename BASE>
std::ostream& trajectory::Adaptive<A,BASE>::displayParameters_v(std::ostream& os) const
{
  return evolved_->displayParameters(os)<<"Trajectory Parameters: epsRel="<<evolved_->getEpsRel()<<" epsAbs="<<evolved_->getEpsAbs()<<std::endl;
}

template<typename A, typename BASE>
cpputils::iarchive& trajectory::Adaptive<A,BASE>::readState_v(cpputils::iarchive& iar)
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

template<typename A, typename BASE>
cpputils::oarchive& trajectory::Adaptive<A,BASE>::writeState_v(cpputils::oarchive& oar) const
{
  meta_.trajectoryID = trajectoryID(); // it is set here rather than @ construction as it is not good to call virtual functions @ construction
  AdaptiveIO<A>::writeState(oar);
  return writeStateMore_v(oar);
}


template<typename A, typename BASE>
void trajectory::Adaptive<A,BASE>::step(double deltaT)
{
  step_v(deltaT);
  if (logLevel_>3)
    this->getLogStreamDuringRun()<<"Number of failed steps in this timestep: "<<evolved_->nFailedStepsLast()<<std::endl;
}



#endif // CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED
