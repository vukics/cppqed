// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Archive.h"
#include "CommentingStream.h"
#include "ODE.h"
#include "Version.h"

#include <iostream>
#include <fstream>
#include <queue>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>



namespace cppqedutils {


template <typename T> concept time_keeper = requires (const T& t) { { getTime(t) } -> std::convertible_to<double>; };

template <typename T> concept adaptive_time_keeper = adaptive_timestep_keeper<T> && time_keeper<T>;


// single adaptive step of at most deltaT
template <typename T>
concept adaptive_steppable = adaptive_time_keeper<T> && requires (T&& t, double deltaT) { { step(t,deltaT) } -> std::convertible_to<LogTree>; };



/// advances for exactly time `deltaT`
LogTree advance(adaptive_steppable auto& traj, double deltaT)
{
  json::array res;
  double endTime=getTime(traj)+deltaT;
  while (double dt=endTime-getTime(traj)) res.push_back(step(traj,dt)) ;
  return {{"advance",res}};
}


namespace trajectory {


// the basic trajectory concept
template <typename T>
concept uniform_step = adaptive_time_keeper<T> &&  logger<T> && requires (T&& t, double deltaT)
  {
    { advance(t,deltaT) } -> std::convertible_to<LogTree>;
    { temporalDataPoint(t) } -> temporal_data_point;
    { dataStreamKey(t) } -> std::convertible_to<LogTree>;
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
LogTree advanceTo(uniform_step auto& traj, double t) { return advance(traj,t-getTime(traj)); }


using StreamSwitch = std::bitset<4>; // stream intro,  stream outro, stream log underway, serialize 1st state


/// Parameters of run() control
/** This could be more finely grained by factoring out autoStopping-control, but the present solution should be sufficient for most puposes */
template <typename BASE = Empty>
struct Pars : BASE
{
  double T; ///< endtime of the run
  unsigned dc;
  double Dt;
  size_t NDt; ///< number of deltaT intervals in \link trajectory::Trajectory::run deltaT-mode\endlink
  std::string ofn, initialFileName;
  
  int precision=6;

  StreamSwitch streamSwitch;

  unsigned sdf;

  double
    autoStopEpsilonRel, ///< relative precision for autostopping
    autoStopEpsilonAbs; ///< absolute precision for autostopping (everything below is not considered)

  unsigned autoStopRepetition; ///< number of streamed lines repeated within relative precision before autostopping – 0 means no autostopping
  
  Pars(popl::OptionParser& op) : BASE{op}
  {
    using ::parameters::_;
    add(op,"Trajectory",  //Title(add(add(add(add(add(add(add(add(add(add(add(op,
     _("T","Simulated time",1.,T),
     _("dc","Number of steps between two streamings",10,dc),
     _("Dt","Timestep between two streamings",.1,Dt),
     _("NDt","Number of steps in Dt mode",0,NDt),
     _("o","Output file name for Trajectory, when empty, cout","",ofn),
     _("initialFileName","Trajectory initial file name","",initialFileName),
//   _("precision","General precision of output",formdouble::Zero(FormDouble::defaultPrecision)),
     _("streamSwitch","stream intro,  stream outro, stream log underway, serialize 1st state",StreamSwitch("1011"),streamSwitch),
     _("sdf","State output frequency",0,sdf),
     _("autoStopEpsilonRel","Relative precision for autostopping",ode::epsRelDefault,autoStopEpsilonRel),
     _("autoStopEpsilonAbs","Absolute precision for autostopping",ode::epsAbsDefault,autoStopEpsilonAbs),
     _("autoStopRepetition","Number of streamed lines repeated within relative precision before autostopping",0,autoStopRepetition)
    );

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
concept observer = uniform_step<TRAJ> && requires (T&& t, const TRAJ& traj) { t(temporalDataPoint(traj)); };

/// Generic implementation of AutostopHandler
/** 
 * Stops if the same value occurs in the output autoStopRepetition times within the specified absolute and relative precisions
 * Assumption is that TDP is a standard-conforming range
 */
template<temporal_data_point TDP>
class AutostopHandlerGeneric
{
public:
  AutostopHandlerGeneric(double autoStopEpsilon, double autoStopEpsilonAbs, unsigned autoStopRepetition) 
    : autoStopEpsilon_{autoStopEpsilon}, autoStopEpsilonAbs_{autoStopEpsilonAbs}, autoStopRepetition_{autoStopRepetition}, averages_{}, queue_{}
    {
      std::ranges::fill(averages_,0.); // this assumes that the value_type of TDP can be converted from double
    }

  void operator()(TDP tdp)
  {
    if (!autoStopRepetition_) return;

    // if there are 0 elements in tdp, replace them with autoStopEpsilonAbs_
    {
      auto replaceZeros = [&] <temporal_data_point T> (T& t, const auto& lambda) {
        if constexpr (std::same_as<T,std::valarray<dcomp>>) t[abs(t)<autoStopEpsilonAbs_]=autoStopEpsilonAbs_;
        else for (auto& v : tdp) lambda(t,lambda);
      };
      replaceZeros(tdp,replaceZeros);
    }

    // This means that the present call is the first
    if (!averages_.size()) {
      averages_.swap(tdp);
      return;
    }
    
    // Queue is fully filled
    if (queue_.size()==autoStopRepetition_) {
      if ( abs(averages_-tdp)/(abs(averages_)+abs(tdp)) < autoStopEpsilon_ ) throw StoppingCriterionReachedException{};
      averages_+=(tdp-queue_.front())/double(autoStopRepetition_); // update the averages set for next step
      queue_.pop(); // remove obsolate first element of the queue
    }
    // queue is not yet fully filled – build the initial averages set to compare against
    else {
      averages_*=double(queue_.size())/double(queue_.size()+1);
      averages_+=tdp/double(queue_.size()+1);
    }

    queue_.push(tdp); // place the new item to the back of the queue

  }

private:
  const double autoStopEpsilon_, autoStopEpsilonAbs_;
  const unsigned autoStopRepetition_;

  TDP averages_;
  std::queue<TDP> queue_;

};


auto observerNoOp = [] (const auto&) {};


template<temporal_data_point TDP>
using DataStream=std::list<std::tuple<double,double,TDP>>;



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


inline std::ostream& stream(double tdp, std::ostream& os) {return os<<tdp;}
inline std::ostream& stream(dcomp  tdp, std::ostream& os) {return os<<tdp;}

inline std::ostream& stream(hana::tuple<> tdp, std::ostream& os) {return os;}

std::ostream& stream(const temporal_data_point auto& tdp, std::ostream& os)
{
  size_t n{0};
  hana::for_each( tdp, [&] (const auto& v) { stream(v,os) << (++n != hana::size(tdp) ? "\t" : "") ; } );
  return os;
}


template <typename T, typename TRAJ>
concept data_streamer = uniform_step<TRAJ> && requires (T&& t, const TRAJ& traj, std::ostream& os) {
  { t(traj,os) } -> std::convertible_to<decltype(temporalDataPoint(traj))> ;
};


const auto dataStreamerDefault = [] (const uniform_step auto& traj, std::ostream& os) {
  auto tdp{temporalDataPoint(traj)};
  stream(tdp, os)<<std::endl; // Note: endl flushes the buffer
  return tdp;
};


/// The most general run function
template < RunLengthType RLT, StreamFreqType SFT, uniform_step TRAJ, data_streamer<TRAJ> TDS = decltype(dataStreamerDefault) >
requires ( ( SFT==StreamFreqType::DT_MODE || adaptive<TRAJ> ) && ( RLT==RunLengthType::T_MODE || SFT==StreamFreqType::DT_MODE ) )
auto
run(TRAJ&& traj, ///< the trajectory to run
    std::conditional_t<bool(RLT),size_t,double> length, ///< length of run
    std::conditional_t<bool(SFT),size_t,double> streamFreq, ///< interval between two \link Trajectory::stream streamings\endlink
    unsigned stateStreamFreq, ///< number of \link Trajectory::stream streamings\endlink between two \link Trajectory::writeState state streamings\endlink
    const std::string& trajectoryFileName, ///< name of the output file for \link Trajectory::stream streams\endlink — if empty, stream to standard output
    const std::string& initialFileName, ///< name of file containing initial condition state for the run
    int precision, ///< governs the overall precision (number of digits) of outputs in \link Trajectory::stream streamings\endlink
    StreamSwitch streamSwitch, bool doStreaming, ///< If false, all trajectory output is redirected to a null-stream
    bool returnStreamedArray, ///< If true, the streamed array is stored and returned by the function
    observer<TRAJ> auto&& observer,
    const TDS& tds = dataStreamerDefault)
{
  auto streamWrapper=[&] (std::ostream& os) {return tds(traj,os<<getTime(traj)<<" "<<getDtDid(traj)<<"\t");};

  using namespace std;

  DataStream<decltype(temporalDataPoint(traj))> res;
  
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

  struct : std::streambuf { int overflow(int c) override { return c; } } nullBuffer;

  const auto outstream{
    doStreaming ? (
      streamToFile ?
      static_pointer_cast<ostream>(make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app)) :
      shared_ptr<ostream>(&cout,[](auto*){}) // since cout is a system-wide object, this should be safe
      ) : make_shared<ostream>(&nullBuffer)
  }; // regulates the deletion policy
  
  if (outstream->fail()) throw runtime_error("Trajectory stream opening error: "+trajectoryFileName);
  
  ostream& os=*outstream;

  CommentingStream logStream{outstream};
  
  ///////////////////////
  // Writing introduction
  ///////////////////////

  if (streamSwitch[0]) {
    if (!continuing) {
      if (parsedCommandLine!="") logStream<<parsedCommandLine<<endl<<endl;
      logStream<<versionHelper()<<endl<<logIntro(traj)<<endl<<"Key to data:\nTrajectory\n 1. time\n 2. dtDid\n"<<dataStreamKey(traj)
        <<endl<<endl<<"Run Trajectory up to time "<<timeToReach
        <<" -- Stream period: "<<streamFreq<< (SFT==StreamFreqType::DT_MODE ? "" : " timestep") <<endl<<endl;
    }
    else
      logStream<<"Continuing from time "<<getTime(traj)<<" up to time "<<timeToReach<<endl;
  }

  if (!timeToReach) {streamWrapper(os); return res;}

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  const std::shared_ptr<ostream> ofs = !streamToFile ? std::make_shared<ofstream>() : openStateFileWriting(stateFileName);

  bool
    stateSaved=false,  // signifies whether the state has already been saved for the actual time instant of the trajectory
    tdpStreamed=false; // signifies whether the temporal data point has already been streamed ”

  try {

    for (long count=0, stateCount=0; (RLT==RunLengthType::T_MODE ? getTime(traj)<length : count<=length); ++count) {

      if (count) {
        // advance trajectory
        LogTree lt;
        if constexpr (SFT==StreamFreqType::DC_MODE) {lt=step(traj,length-getTime(traj));}
        else {
          if constexpr (RLT==RunLengthType::T_MODE) lt=advance(traj,std::min(streamFreq,length-getTime(traj)));
          else lt=advance(traj,streamFreq);
        }
        if (streamSwitch[2] && size(lt)) logStream<<lt<<endl;
        stateSaved=tdpStreamed=false;
      }

      if (!count || [&]() {
        if constexpr (SFT==StreamFreqType::DT_MODE) return true;
        else return !(count%streamFreq); // here, we still use a lambda because this doesn’t compile if streamFreq is double
      }() ) {

        if (
            stateStreamFreq && 
            !(stateCount%stateStreamFreq) && 
            (stateCount || (streamSwitch[3] && !continuing))
           )
        {
          writeViaSStream(traj,ofs);
          stateSaved=true;
        }
        ++stateCount;

        if (count || !continuing) {
          tdpStreamed=true;
          auto streamReturn{streamWrapper(os)};
          if (returnStreamedArray) res.emplace_back(getTime(traj),getDtDid(traj),streamReturn);
          observer(streamReturn);
        }
      }
    } // end of main for loop

  } catch (const StoppingCriterionReachedException& except) {logStream<<"Stopping criterion has been reached"<<endl;}
  
  if (!tdpStreamed) streamWrapper(os); // Stream at the end instant if stream has not happened yet

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  if (streamSwitch[1]) logStream<<"OUTRO: "<<logOutro(traj)<<endl;
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
template<typename AH, trajectory::uniform_step TRAJ, typename PB> requires trajectory::observer<AH,TRAJ>
auto run(TRAJ&& traj, const trajectory::Pars<PB>& p, AH&& ah, bool doStreaming=true, bool returnStreamedArray=false)
{
  using namespace trajectory;
  if constexpr (adaptive<TRAJ>) {
    if (p.dc) return trajectory::run<RunLengthType::T_MODE,StreamFreqType::DC_MODE>(
      traj,p.T,p.dc,p.sdf,p.ofn,p.initialFileName,p.precision,p.streamSwitch,doStreaming,returnStreamedArray,ah);
    else if (!p.Dt) throw std::runtime_error("Nonzero dc or Dt required in trajectory::run");
  }
  if (!p.Dt) throw std::runtime_error("Nonzero Dt required in trajectory::run");
  if (p.NDt) return trajectory::run<RunLengthType::NDT_MODE,StreamFreqType::DT_MODE>(
    traj,p.NDt,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,p.streamSwitch,doStreaming,returnStreamedArray,ah);
  else return trajectory::run<RunLengthType::T_MODE,StreamFreqType::DT_MODE>(
    traj,p.T,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,p.streamSwitch,doStreaming,returnStreamedArray,ah);
}


template<trajectory::uniform_step TRAJ, typename PB>
auto run(TRAJ&& traj, const trajectory::Pars<PB>& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return run(traj,p,trajectory::AutostopHandlerGeneric<decltype(temporalDataPoint(traj))>(p.autoStopEpsilonRel,p.autoStopEpsilonAbs,p.autoStopRepetition),
             doStreaming,returnStreamedArray);
}


} // cppqedutils

