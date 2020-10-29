// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_TRAJECTORY_TCC_INCLUDED

#include "Trajectory.h"

#include "CommentingStream.h"
#include "Evolved.tcc"
#include "FormDouble.h"
#include "IO_Manip.h"
#include "Version.h"

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/stream.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>



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

template <typename DA> inline bool doContinue(const Trajectory<DA>& traj, double time, long      ) {return traj.getTime()<time;}
template <typename DA> inline bool doContinue(const Trajectory<DA>&     , long length, long count) {return count<=length      ;}

template <typename DA> inline void advance(Trajectory<DA>& traj, long       , double deltaT) {traj.evolve(deltaT);}
template <typename DA> inline void advance(Trajectory<DA>& traj, double time, double deltaT) {traj.evolve(std::min(deltaT,time-traj.getTime()));}

template<typename A, typename BASE>
inline void advance(Adaptive<A,BASE>& traj, double time, int) {traj.step(time-traj.getTime());}

inline bool doDisplay(long      , double         ) {return true;}
inline bool doDisplay(long count, int displayFreq) {return !(count%displayFreq);}

inline auto writeTimestep(int   ) {return " timestep";}
inline auto writeTimestep(double) {return ""         ;}

inline double endTime(long    nDt, double dt, double currentTime) {return nDt*dt+currentTime;}
inline double endTime(double time, double   , double            ) {return time              ;}

} // runTraits


template<typename DA>
bool restoreState(Trajectory<DA>& traj, const std::string& trajectoryFileName, const std::string& stateFileName, const std::string& initialFileName)
{
  if (trajectoryFileName!="") {

    std::ifstream trajectoryFile(trajectoryFileName.c_str());

    if (trajectoryFile.is_open() && (trajectoryFile.peek(), !trajectoryFile.eof()) ) {
      auto stateFile{openStateFileReading(stateFileName)};
      while ( (stateFile->peek(), !stateFile->eof()) ) readViaSStream(traj,stateFile.get());
      return true;
    }

  }

  if (initialFileName!="") {
    auto initialFile{openStateFileReading(initialFileName)};
    while ( (initialFile->peek(), !initialFile->eof()) ) readViaSStream(traj,initialFile.get());
  }

  return false;
}

} } // trajectory::details


template<typename T, typename L, typename D, typename AutostopHandler>
trajectory::TrajectoryDisplayList<typename T::DisplayedArray> 
trajectory::details::run(T& traj, L length, D displayFreq, unsigned stateDisplayFreq,
                         const std::string& trajectoryFileName, const std::string& initialFileName,
                         int precision, bool displayInfo, bool firstStateDisplay,
                         const std::string& parsedCommandLine,
                         bool doStreaming, bool saveDisplayedArray,
                         AutostopHandler&& autostopHandler
                        )
{
  using namespace std; using namespace runTraits; using namespace cpputils;

  TrajectoryDisplayList<typename T::DisplayedArray> res;
  
  ////////////////////////////////////////////////
  // Determining i/o streams, eventual state input
  ////////////////////////////////////////////////

  static const string stateExtension(".state");
  const string stateFileName(trajectoryFileName+stateExtension);
  
  const bool
    streamToFile=(trajectoryFileName!=""),  
    continuing=restoreState(traj,trajectoryFileName,stateFileName,initialFileName);
  
  const double timeToReach=endTime(length,displayFreq,traj.getTime());
    
  if (timeToReach && timeToReach<=traj.getTime()) return res;

  const std::shared_ptr<ostream> outstream(!streamToFile ?
                                           ( doStreaming ? 
                                             std::shared_ptr<ostream>(&cout,[](auto*){}) : // since cout is a system-wide object, this should be safe
                                             static_pointer_cast<ostream>(std::make_shared<boost::iostreams::stream<boost::iostreams::null_sink>>(boost::iostreams::null_sink{}) ) ) : 
                                           static_pointer_cast<ostream>(std::make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app))); // regulates the deletion policy
  
  if (outstream->fail()) throw std::runtime_error("Trajectory file opening error: "+trajectoryFileName);
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

  if (!timeToReach) {traj.display(os,precision); return res;}

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  const std::shared_ptr<ostream> ofs = !streamToFile ? std::make_shared<ofstream>() : openStateFileWriting(stateFileName);

  bool
    stateSaved=false,   // signifies whether the state has already been saved for the actual time instant of the trajectory
    evsDisplayed=false; // signifies whether the expectation values have already been displayed ”

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

        if (count || !continuing) {
          evsDisplayed=true;
          auto displayReturn{traj.display(os,precision)};
          if (saveDisplayedArray) res.emplace_back(traj.getTime(),traj.getDtDid(),get<1>(displayReturn));
          autostopHandler(get<1>(displayReturn));
        }
      }
    }

  } catch (const StoppingCriterionReachedException& except) {CommentingStream(os)<<"Stopping criterion has been reached"<<endl;}
  
  if (!evsDisplayed) traj.display(os,precision); // Display at the end instant if display has not happened yet

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  {CommentingStream temp(os); temp<<setprecision(formdouble::actualPrecision(precision)); traj.logOnEnd(temp);}
  if (!stateSaved) writeViaSStream(traj,ofs.get());
  
  return res;
  
}



template<typename A>
trajectory::AdaptiveIO<A>::AdaptiveIO(typename EvolvedIO::Ptr evolvedIO)
  : meta_(cpputils::TypeID<A>::value,
          SerializationMetadata::ARRAY_ONLY,
          cpputils::Rank<A>::value),
    evolvedIO_(evolvedIO)
{}

template<typename A>
cpputils::iarchive& trajectory::AdaptiveIO<A>::readState(cpputils::iarchive& iar)
{
  bool dimension_check = meta_.trajectoryID != SerializationMetadata::ARRAY_ONLY;
  iar & meta_;
  if (meta_.rank!=cpputils::Rank<A>::value)
    throw std::runtime_error("Rank mismatch in trajectory::AdaptiveIO::readState");
  std::vector<size_t> dims = cpputils::dimensions(evolvedIO_->getA());
  iar & *evolvedIO_;
  if (dimension_check && dims != cpputils::dimensions(evolvedIO_->getA()))
      throw std::runtime_error("Dimensions mismatch in trajectory::AdaptiveIO::readState");
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
      throw std::runtime_error("Trajectory mismatch in trajectory::Adaptive::readState_v");
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
