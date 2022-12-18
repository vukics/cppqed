// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Archive.h"
#include "CommentingStream.h"
#include "FormDouble.h"
#include "ODE.h"
#include "Version.h"

#include <boost/hana.hpp>

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/stream.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>



namespace hana=boost::hana;


namespace cppqedutils {


template <typename T> concept time_keeper = requires (const T& t) { { getTime(t) } -> std::convertible_to<double>; };

template <typename T> concept adaptive_time_keeper = adaptive_timestep_keeper<T> && time_keeper<T>;


// single adaptive step of at most deltaT
template <typename T>
concept adaptive_steppable = adaptive_time_keeper<T> && requires (T&& t, double deltaT, std::ostream& os) { step(t,deltaT,os); };



/// advances for exactly time `deltaT`
template<adaptive_steppable T>
void advance(T& traj, double deltaT, std::ostream& logStream)
{
  double endTime=getTime(traj)+deltaT;
  while (double dt=endTime-getTime(traj)) step(traj,dt,logStream);
}


namespace trajectory {
  

// the basic trajectory concept
template <typename T>
concept uniform_step =
  adaptive_time_keeper<T> && 
  intro_outro_streamer<T> && 
  requires (T&& t, double deltaT, std::ostream& os, int precision)
  {
    advance(t,deltaT,os);
    { streamKey(t,os) } -> std::convertible_to<std::ostream&>;
    { stream(t,os,precision) } -> std::convertible_to<typename std::decay_t<T>::StreamedArray>;
  };

template <typename T>
concept adaptive = uniform_step<T> && adaptive_steppable<T>;


template <typename T, typename Archive>
concept archiving = requires (T&& t, Archive& ar) {
  { stateIO(t,ar) } -> std::convertible_to<Archive&>;
  // { readFromArrayOnlyArchive(t,ar) } -> std::convertible_to<Archive&>;
  // this latter only holds for iarchives, so we should split the concept into input/output archiving if we want to include this
};


enum struct StreamFreqType {DT_MODE=0, DC_MODE=1};

enum struct RunLengthType {T_MODE=0, NDT_MODE=1};


/// advances up to exactly time `t` \copydetails advance
template<uniform_step Trajectory>
void advanceTo(Trajectory& traj, double t, std::ostream& logStream) { advance(traj,t-getTime(traj),logStream); }



/// Parameters corresponding to the different versions of run()
template <typename BASE = Empty>
struct Pars : BASE
{
  double T; ///< endtime of the run
  unsigned dc;
  double Dt;
  size_t NDt; ///< number of deltaT intervals in \link trajectory::Trajectory::run deltaT-mode\endlink
  std::string ofn, initialFileName;
  
  int precision=6;

  bool streamInfo=true, firstStateStream=true;

  unsigned sdf;

  double
    autoStopEpsilon, ///< relative precision for autostopping
    autoStopEpsilonAbs; ///< absolute precision for autostopping (everything below is not considered)

  unsigned autoStopRepetition; ///< number of streamed lines repeated within relative precision before autostopping – 0 means no autostopping
  
  Pars(popl::OptionParser& op) : BASE{op}
  {
    addTitle(add(add(add(add(add(add(op,
     "initialFileName","Trajectory initial file name",std::string{},&initialFileName),
     "o","Output file name for Trajectory, when empty, cout",std::string{},&ofn),
     "NDt","Number of steps in Dt mode",size_t{0},&NDt),
     "Dt","Timestep between two streamings",.1,&Dt),
     "dc","Number of steps between two streamings",10u,&dc),
     "T","Simulated time",1.,&T),
     "Trajectory");
/*    
    streamInfo(p.add("streamInfo",mod,"Whether to stream header for trajectories",true)),
    firstStateStream(p.add("firstStateStream",mod,"Streams trajectory state at startup",true)),
    sdf(p.add("sdf",mod,"State output frequency",0u)),
    autoStopEpsilon(p.add("autoStopEpsilon",mod,"Relative precision for autostopping",ode_engine::epsRelDefault)),
    autoStopEpsilonAbs(p.add("autoStopEpsilonAbs",mod,"Absolute precision for autostopping",ode_engine::epsAbsDefault)),
    autoStopRepetition(p.add("autoStopRepetition",mod,"Number of streamed lines repeated within relative precision before autostopping",0u)), */
  }

};


/// A heuristic determination of the inital timestep from the highest frequency of a physical system.
inline double initialTimeStep(double highestFrequency) {return 1./(10.*highestFrequency);}
  

/// Type-erased aggregate of information about a trajectory-state archive
/**
 * \par Rationale
 * Each archive is self-contained, so that it contains its own metadata.
 * Then, the state files can eventually even be split into several archives, each containing a single self-contained state.
 */
struct SerializationMetadata
{ 
  explicit SerializationMetadata(std::string type=UNSPECIFIED, std::string id=UNSPECIFIED, int r=0)
    : rank(r), typeID(type), trajectoryID(id) {};
  
  bool operator==(const SerializationMetadata&) const = default;
  
  int rank=0;
  std::string typeID=UNSPECIFIED, trajectoryID=UNSPECIFIED;

  int protocolVersion=0;  

  inline static const std::string UNSPECIFIED = "Unspecified" ;
  inline static const std::string ARRAY_ONLY = "ArrayOnly" ;

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & protocolVersion & rank & typeID & trajectoryID;}

};


template <uniform_step Trajectory>
struct MakeSerializationMetadata;


struct StoppingCriterionReachedException {};



template <typename T, typename TRAJ>
concept autostop_handler = uniform_step<TRAJ> && requires (T&& t, const typename std::decay_t<TRAJ>::StreamedArray& sa) { t(sa); };

/// Generic implementation of AutostopHandler
/**
 * Assumes quite a lot about the implicit interface of SA, this could be relieved with indirections
 * TODO: should rely on boost::odeint state algebra
 */
template<typename SA>
class AutostopHandlerGeneric
{
private:
  typedef SA Averages;
  
public:
  AutostopHandlerGeneric(double autoStopEpsilon, double autoStopEpsilonAbs, unsigned autoStopRepetition) 
    : autoStopEpsilon_{autoStopEpsilon}, autoStopEpsilonAbs_{autoStopEpsilonAbs}, autoStopRepetition_{autoStopRepetition}, averages_{}, queue_{} {}

  void operator()(const SA& streamedArray)
  {
    SA buffer{streamedArray};
    
    for (auto& v : buffer) if (std::abs(v)<autoStopEpsilonAbs_) v=autoStopEpsilonAbs_;
    
    if (!autoStopRepetition_) return;
    
    if (!averages_.size()) { // This means that no stream has yet occured: the size of averages_ must be determined
      averages_.resize(buffer.size());
      averages_=0.;
    }
    else {
      if (queue_.size()==autoStopRepetition_) {
        if (max(abs(averages_-buffer)/(abs(averages_)+abs(buffer)))<autoStopEpsilon_) throw StoppingCriterionReachedException();
        averages_=averages_+(buffer-queue_.front())/double(autoStopRepetition_); // update the averages set for next step
        queue_.pop(); // remove obsolate first element of the queue
      }
      else averages_=(double(queue_.size())*averages_+buffer)/double(queue_.size()+1); // build the initial averages set to compare against

      queue_.push(buffer); // place the new item to the back of the queue

    }

  }

private:
  const double autoStopEpsilon_, autoStopEpsilonAbs_;
  const unsigned autoStopRepetition_;

  SA averages_;
  std::queue<SA> queue_;

};


const struct
{
  template <typename SA>
  void operator()(const SA&) const {}
} autostopHandlerNoOp;


template<typename SA>
using TemporalStreamedArray=std::list<std::tuple<double,double,SA>>;



/// Opens a file with the given filename for reading and returns an associated std::istream object.
/**
 * If C++QED is compiled with compression support (boost iostreams library needed), then
 * the statefile can be in bzip2-format, the returned istream will handle this automatically.
 *
 * If the file cannot be opened, an exception is raised.
 */
std::shared_ptr<std::istream> openStateFileReading(const std::string &filename);


/// Opens a file with the given filename for writing (appending) and returns an associated std::ostream object.
/**
 * If C++QED is compiled with compression support (boost iostreams library needed), then
 * the statefile can be in bzip2-format, the returned ostream will handle this automatically.
 * If the file does not yet exist, bzip2 is used if supported.
 *
 * If the file cannot be opened, an exception is raised.
 */
std::shared_ptr<std::ostream> openStateFileWriting(const std::string &filename, const std::ios_base::openmode mode=std::ios_base::app);


/// TODO: Read/WriteState and read/writeViaSStream are perhaps not necessary anymore if the interface to numpy is done with JSON
template <archiving<iarchive> Trajectory>
struct ReadState
{
  static iarchive& _(Trajectory& traj, iarchive& iar)
  {
    SerializationMetadata sm;
    iar & sm;
    if (sm.trajectoryID==SerializationMetadata::ARRAY_ONLY) return readFromArrayOnlyArchive(traj,iar);
    else {
      if (sm != MakeSerializationMetadata<Trajectory>::_()) throw std::runtime_error("Trajectory mismatch in cppqedutils::trajectory::ReadState");
      return stateIO(traj,iar);
    }
  }
};


template <archiving<oarchive> Trajectory>
struct WriteState
{
  static oarchive& _(Trajectory& traj, // cannot be const, because traj.stateIO is non-const
                     oarchive& oar)
  {
    return stateIO(traj, oar & MakeSerializationMetadata<Trajectory>::_());
  }  
};


template<archiving<iarchive> Trajectory>
void readViaSStream(Trajectory& traj, std::shared_ptr<std::istream> ifs)
{
  using namespace std;
  istringstream iss(ios_base::binary);
  {
    string buffer;
    streamsize n; *ifs>>n; buffer.resize(n);
    ifs->read(&buffer[0],n);
    iss.str(buffer);
  }
  iarchive stateArchive(iss);
  ReadState<Trajectory>::_(traj,stateArchive);
}


template<archiving<oarchive> Trajectory>
void writeViaSStream(Trajectory& traj, // cannot be const, because traj.stateIO is non-const
                     std::shared_ptr<std::ostream> ofs)
{
  using namespace std;
  if (ofs) {
    ostringstream oss(ios_base::binary);
    oarchive stateArchive(oss);
    WriteState<Trajectory>::_(traj,stateArchive);
    {
      const string& buffer=oss.str();
      *ofs<<buffer.size(); ofs->write(&buffer[0],buffer.size());
    }
  }
}



/// The most general run function
template<RunLengthType RLT, StreamFreqType SFT, uniform_step TRAJ, autostop_handler<TRAJ> AH >
requires ( ( SFT==StreamFreqType::DT_MODE || adaptive<TRAJ> ) && ( RLT==RunLengthType::T_MODE || SFT==StreamFreqType::DT_MODE ) )
TemporalStreamedArray<typename std::decay_t<TRAJ>::StreamedArray>
run(TRAJ&& traj, ///< the trajectory to run
    std::conditional_t<bool(RLT),size_t,double> length, ///< length of run
    std::conditional_t<bool(SFT),size_t,double> streamFreq, ///< interval between two \link Trajectory::stream streamings\endlink
    unsigned stateStreamFreq, ///< number of \link Trajectory::stream streamings\endlink between two \link Trajectory::writeState state streamings\endlink
    const std::string& trajectoryFileName, ///< name of the output file for \link Trajectory::stream streams\endlink — if empty, stream to standard output
    const std::string& initialFileName, ///< name of file containing initial condition state for the run
    int precision, ///< governs the overall precision (number of digits) of outputs in \link Trajectory::stream streamings\endlink
    bool streamInfo, //< governs whether a \link Trajectory::streamIntro header\endlink is streamed at the top of the output
    bool firstStateStream, ///< governs whether the state is streamed at time zero (important if \link Trajectory::writeState state streaming\endlink is costly)
    bool doStreaming, ///< If false, all trajectory output is redirected to a null-stream
    bool returnStreamedArray, ///< If true, the streamed array is stored and returned by the function
    AH&& autostopHandler
   )
{
  auto streamWrapper=[&] (std::ostream& os)
  {
    const FormDouble fd(formdouble::positive(precision));
    auto res{stream(traj, os<<fd(getTime(traj))<<fd(getDtDid(traj)) , precision)};
    os<<std::endl; // Note: endl flushes the buffer
    return res;
  };
 
  using namespace std;

  TemporalStreamedArray<typename decay_t<TRAJ>::StreamedArray> res;
  
  ////////////////////////////////////////////////
  // Determining i/o streams, eventual state input
  ////////////////////////////////////////////////

  const string stateFileName{trajectoryFileName+".state"};
  
  const bool
    streamToFile=(trajectoryFileName!=""),
    continuing=[&]() {
      if (trajectoryFileName!="") {
        ifstream trajectoryFile{trajectoryFileName.c_str()};
        if (trajectoryFile.is_open() && (trajectoryFile.peek(), !trajectoryFile.eof()) ) {
          auto stateFile{openStateFileReading(stateFileName)};
          while ( (stateFile->peek(), !stateFile->eof()) ) readViaSStream(traj,stateFile);
          return true;
        }
      }
      if (initialFileName!="") {
        auto initialFile{openStateFileReading(initialFileName)};
        while ( (initialFile->peek(), !initialFile->eof()) ) readViaSStream(traj,initialFile);
      }
      return false;
    }();
    
  const double timeToReach = (RLT==RunLengthType::T_MODE ? length : getTime(traj)+length*streamFreq);
  
  if (timeToReach && timeToReach<=getTime(traj)) return res;

  const shared_ptr<ostream> outstream{
    !streamToFile ?
    ( doStreaming ? 
      shared_ptr<ostream>(&cout,[](auto*){}) : // since cout is a system-wide object, this should be safe
      static_pointer_cast<ostream>(make_shared<boost::iostreams::stream<boost::iostreams::null_sink>>(boost::iostreams::null_sink{}) ) ) : 
    static_pointer_cast<ostream>(make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app))}; // regulates the deletion policy
  
  if (outstream->fail()) throw runtime_error("Trajectory stream opening error: "+trajectoryFileName);
  
  CommentingStream logStream{outstream};
  
  ostream& os=*outstream;
  os<<setprecision(formdouble::actualPrecision(precision));

  CommentingStream commentingStream{outstream}; commentingStream<<setprecision(formdouble::actualPrecision(precision));
  
  ///////////////////////
  // Writing introduction
  ///////////////////////

  if (streamInfo) {
    if (!continuing) {
      if (parsedCommandLine!="") commentingStream<<parsedCommandLine<<endl<<endl;
      streamKey(traj,streamIntro(traj,commentingStream<<versionHelper())<<std::endl<<"Key to data:\nTrajectory\n 1. time\n 2. dtDid\n")
        <<endl<<"Run Trajectory up to time "<<timeToReach
        <<" -- Stream period: "<<streamFreq<< (SFT==StreamFreqType::DT_MODE ? "" : " timestep") <<endl<<endl;
    }
    else
      commentingStream<<"Continuing from time "<<getTime(traj)<<" up to time "<<timeToReach<<endl;
  }

  if (!timeToReach) {streamWrapper(os); return res;}

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  const std::shared_ptr<ostream> ofs = !streamToFile ? std::make_shared<ofstream>() : openStateFileWriting(stateFileName);

  bool
    stateSaved=false,  // signifies whether the state has already been saved for the actual time instant of the trajectory
    arrayStreamed=false; // signifies whether the expectation values have already been streamed ”

  try {

    for (long count=0, stateCount=0; (RLT==RunLengthType::T_MODE ? getTime(traj)<length : count<=length); ++count) {

      if (count) {
        // advance trajectory
        if constexpr (SFT==StreamFreqType::DC_MODE) {step(traj,length-getTime(traj),logStream);}
        else {
          if constexpr (RLT==RunLengthType::T_MODE) advance(traj,std::min(streamFreq,length-getTime(traj)),logStream);
          else advance(traj,streamFreq,logStream);
        }
        stateSaved=arrayStreamed=false;
      }

      if (!count || [&]() {
        if constexpr (SFT==StreamFreqType::DT_MODE) return true;
        else return !(count%streamFreq); // here, we still use a lambda because this doesn’t compile if streamFreq is double
      }() ) {

        if (
            stateStreamFreq && 
            !(stateCount%stateStreamFreq) && 
            (stateCount || (firstStateStream && !continuing))
           )
        {
          writeViaSStream(traj,ofs);
          stateSaved=true;
        }
        ++stateCount;

        if (count || !continuing) {
          arrayStreamed=true;
          auto streamReturn{streamWrapper(os)};
          if (returnStreamedArray) res.emplace_back(getTime(traj),getDtDid(traj),streamReturn);
          autostopHandler(streamReturn);
        }
      }
    } // end of main for loop

  } catch (const StoppingCriterionReachedException& except) {commentingStream<<"Stopping criterion has been reached"<<endl;}
  
  if (!arrayStreamed) streamWrapper(os); // Stream at the end instant if stream has not happened yet

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  streamOutro(traj,commentingStream);
  if (!stateSaved) writeViaSStream(traj,ofs);
  
  return res;
  
}


} // trajectory


/// Dispatcher \related Adaptive
/**
 * Since in addition to the Trajectory interface, Adaptive has also the capability to be propagated over a \link Adaptive::step single adaptive timestep\endlink,
 * it is possible to count the individual ODE steps. Running in dc-mode means that a fixed number of adaptive steps are performed between each streaming.
 * Hence, we get denser streamings in time when the timestep is small, that is, when the important things are happening in the dynamics.
 *
 * NDt-mode: runs for a certain number of time intervals deltaT \related Trajectory
 *
 * <b>Rationale:</b> This version of `run` exists to avoid the eventual tiny timestep at the end of the run that might occur with
 * \link Trajectory::run the above version\endlink.
 * This is because e.g. in the case of `deltaT=0.1`, `time=1`, adding up `0.1` ten times numerically does not result in exactly `1`.
 * 
 * For a demonstration, compare the output of
 *
 *     examples/HarmonicOscillatorComplex --dc 0 --Dt 0.1 --T 1
 *
 * with
 *
 *     examples/HarmonicOscillatorComplex --dc 0 --Dt 0.1 --NDt 10
 *
 *
 * - Runs run(Adaptive&, double, int, unsigned, const std::string&, const std::string&, int, bool, bool) if `p.dc` is nonzero (dc-mode)
 * - delegates to run(Trajectory&, const ParsRun& p) otherwise (deltaT-mode)
 * 
 * \note This means that ParsRun::NDt takes precedence over ParsRun::T and ParsRun::dc takes precedence over ParsRun::Dt
 * 
 */
template<typename AH, trajectory::uniform_step TRAJ, typename PB> requires trajectory::autostop_handler<AH,TRAJ>
auto run(TRAJ&& traj, const trajectory::Pars<PB>& p, AH&& ah, bool doStreaming=true, bool returnStreamedArray=false)
{
  using namespace trajectory;
  if constexpr (adaptive<TRAJ>) {
    if (p.dc) return trajectory::run<RunLengthType::T_MODE,StreamFreqType::DC_MODE>(
      std::forward<TRAJ>(traj),p.T,p.dc,p.sdf,p.ofn,p.initialFileName,p.precision,
      p.streamInfo,p.firstStateStream,
      doStreaming,returnStreamedArray,std::forward<AH>(ah));
    else if (!p.Dt) throw std::runtime_error("Nonzero dc or Dt required in trajectory::run");
  }
  if (!p.Dt) throw std::runtime_error("Nonzero Dt required in trajectory::run");
  if (p.NDt) return trajectory::run<RunLengthType::NDT_MODE,StreamFreqType::DT_MODE>(
    std::forward<TRAJ>(traj),p.NDt,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,
    p.streamInfo,p.firstStateStream,
    doStreaming,returnStreamedArray,std::forward<AH>(ah));
  else return trajectory::run<RunLengthType::T_MODE,StreamFreqType::DT_MODE>(
    std::forward<TRAJ>(traj),p.T,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,
    p.streamInfo,p.firstStateStream,
    doStreaming,returnStreamedArray,std::forward<AH>(ah));
}


template<trajectory::uniform_step TRAJ, typename PB>
auto
run(TRAJ&& traj, const trajectory::Pars<PB>& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return run(traj,p,trajectory::AutostopHandlerGeneric<typename std::decay_t<TRAJ>::StreamedArray>(p.autoStopEpsilon,p.autoStopEpsilonAbs,p.autoStopRepetition),
             doStreaming,returnStreamedArray);
}


} // cppqedutils

