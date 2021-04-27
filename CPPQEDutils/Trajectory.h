// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the basic classes of the trajectory-bundle}
#ifndef CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED

#include "Archive.h"
#include "CommentingStream.h"
#include "FormDouble.h"
#include "IO_Manip.h"
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


/// The trajectory-bundle
namespace cppqedutils {


inline auto has_step = hana::is_valid([](auto&& obj) -> decltype(obj.step(1.0,std::clog)) { });
inline auto has_advance = hana::is_valid([](auto&& obj) -> decltype(obj.advance(1.0,std::clog)) { });


/// \name Generic evolution functions
//@{
/// advances for exactly time `deltaT`
/** \tparam Trajectory type of the trajectory to advance. Enabled only if it has `step(double, std::ostream&)` */
template<typename Trajectory>
void advance(Trajectory& traj,
            std::enable_if_t<decltype(has_step(traj))::value,double> deltaT,
            std::ostream& logStream=std::clog)
{
  double endTime=traj.getTime()+deltaT;
  while (double dt=endTime-traj.getTime()) traj.step(dt,logStream);
}


/// advances for exactly time `deltaT`
/** \tparam Trajectory type of the trajectory to advance. Enabled only if it has `advance(double, std::ostream&)` */
template<typename Trajectory>
void advance(Trajectory& traj,
            std::enable_if_t<!decltype(has_step(traj))::value && decltype(has_advance(traj))::value,double> deltaT,
            std::ostream& logStream=std::clog)
{
  traj.advance(deltaT,logStream);
}


/// advances up to exactly time `t` \copydetails advance
template<typename Trajectory>
void advanceTo(Trajectory& traj, double t, std::ostream& logStream=std::clog)
{
  advance(traj,t-traj.getTime(),logStream);
}
//@}

  
namespace trajectory {
  

/// Parameters corresponding to the different versions of run()
template <typename BASE = ode_engine::Pars<>>
struct Pars : BASE
{
  double &T; ///< endtime of the run
  int &dc;
  double &Dt;
  long &NDt; ///< number of deltaT intervals in \link trajectory::Trajectory::run deltaT-mode\endlink
  std::string &ofn, &initialFileName;

  formdouble::Zero &precision; ///< the overall precision of trajectory stream \see FormDouble::overallPrecision

  bool &streamInfo, &firstStateStream;

  unsigned &sdf;

  double
    &autoStopEpsilon, ///< relative precision for autostopping
    &autoStopEpsilonAbs; ///< absolute precision for autostopping (everything below is not considered)

  unsigned &autoStopRepetition; ///< number of streamed lines repeated within relative precision before autostopping – 0 means no autostopping
  
  Pars(parameters::Table& p, const std::string& mod="")
  : BASE{p,mod},
    T(p.addTitle("Trajectory",mod).add("T",mod,"Simulated time",1.)),
    dc(p.add("dc",mod,"Number of steps between two streamings",10)),
    Dt(p.add("Dt",mod,"Timestep between two streamings",.1)),
    NDt(p.add("NDt",mod,"Number of steps in Dt mode",0L)),
    ofn(p.add<std::string>("o",mod,"Output file name for Trajectory, when empty, cout","")),
    initialFileName(p.add<std::string>("initialFileName",mod,"Trajectory initial file name","")),
    precision(p.add("precision",mod,"General precision of output",formdouble::Zero(FormDouble::defaultPrecision))),
    streamInfo(p.add("streamInfo",mod,"Whether to stream header for trajectories",true)),
    firstStateStream(p.add("firstStateStream",mod,"Streams trajectory state at startup",true)),
    sdf(p.add("sdf",mod,"State output frequency",0u)),
    autoStopEpsilon(p.add("autoStopEpsilon",mod,"Relative precision for autostopping",ode_engine::epsRelDefault)),
    autoStopEpsilonAbs(p.add("autoStopEpsilonAbs",mod,"Absolute precision for autostopping",ode_engine::epsAbsDefault)),
    autoStopRepetition(p.add("autoStopRepetition",mod,"Number of streamed lines repeated within relative precision before autostopping",0u)),
    parsedCommandLine_(p.getParsedCommandLine())
    {};

  /// Corresponds to parameters::Table::getParsedCommandLine
  const std::string getParsedCommandLine() const {return *parsedCommandLine_;}

private:
  const parameters::Table::ParsedCommandLine parsedCommandLine_;

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
  
  bool operator==(const SerializationMetadata& v) const {return rank==v.rank && typeID==v.typeID && trajectoryID==v.trajectoryID;}
  bool operator!=(const SerializationMetadata& v) const {return !(*this==v);}
  
  int rank=0;
  std::string typeID=UNSPECIFIED;
  std::string trajectoryID=UNSPECIFIED;

  int protocolVersion=0;  

  inline static const std::string UNSPECIFIED = "Unspecified" ;
  inline static const std::string ARRAY_ONLY = "ArrayOnly" ;

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & protocolVersion & rank & typeID & trajectoryID;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

};


template <typename Trajectory>
struct MakeSerializationMetadata;


struct StoppingCriterionReachedException {};


} // trajectory


namespace trajectory {


/// Generic implementation of AutostopHandler
/**
 * Assumes quite a lot about the implicit interface of SA, this could be relieved with indirections
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
    SA buffer{streamedArray.copy()};
    
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


/// Running in deltaT mode (streams in equal time intervals) for a certain time \related Trajectory
/**
 * This function manifests all the basic features of Adaptive and the whole idea behind the trajectory bundle.
 *
 * A Trajectory can
 * - be \link Trajectory::advance advanced\endlink (propagated in time by given time intervals)
 * - \link Trajectory::readState perform i/o of its entire state\endlink, that is, a bunch of information necessary for resuming a Trajectory from a certain time instant
 * - \link Trajectory::stream stream relevant physical and numerical information\endlink about its actual state at any time (e.g. a set of quantum averages in the case of a quantum trajectory)
 *
 * \note While the entire state can be huge (e.g. the state vector or density operator in the case of a quantum trajectory) the relevant information in an actual numerical experiment
 * is usually much less (a set of quantum averages that doesn’t entirely define the state).
 *
 * Furthermore, a Trajectory can
 * - provide information about its \link Trajectory::getTime time\endlink and \link Trajectory::getDtDid last performed timestep\endlink
 * - \link Trajectory::streamParameters print a header\endlink summarizing its physical and numerical parameters together with a key to the set of relevant physical information streamed during the run
 * - \link Trajectory::logOnEnd print a log\endlink at the end summarizing overall (e.g. time-averaged) physical and numerical data during the run
 *
 * \see Simulated for a full generic implementation of Trajectory together with a small tutorial
 * 
 * \link Trajectory::writeState state streams\endlink into file named `ofn.state`
 *
 */
template<typename TRAJ, // type of trajectory
         typename LENGTH, // type specifying the length of the run
         typename DELTA, // type specifying the frequency of stream
         typename AutostopHandler // should support operator()(const typename T::TemporalStreamedArray &)
         >
TemporalStreamedArray<typename std::decay_t<TRAJ>::StreamedArray>
run(TRAJ&& traj, ///< the trajectory to run
    LENGTH length, ///< length of run
    DELTA streamFreq, ///< interval between two \link Trajectory::stream streamings\endlink
    unsigned stateStreamFreq, ///< number of \link Trajectory::stream streamings\endlink between two \link Trajectory::writeState state streamings\endlink
    const std::string& trajectoryFileName, ///< name of the output file for \link Trajectory::stream streams\endlink — if empty, stream to standard output
    const std::string& initialFileName, ///< name of file containing initial condition state for the run
    int precision, ///< governs the overall precision (number of digits) of outputs in \link Trajectory::stream streamings\endlink
    bool streamInfo, //< governs whether a \link Trajectory::streamParameters header\endlink is streamed at the top of the output
    bool firstStateStream, ///< governs whether the state is streamed at time zero (important if \link Trajectory::writeState state streaming\endlink is costly)
    const std::string& parsedCommandLine, 
    bool doStreaming, ///< If false, all trajectory output is redirected to a null-stream
    bool returnStreamedArray, ///< If true, the streamed array is stored and returned by the function
    AutostopHandler&& autostopHandler
   );


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
template<typename AutostopHandler, typename TRAJ, typename ParsBase>
auto
run(TRAJ&& traj, const trajectory::Pars<ParsBase>& p, AutostopHandler&& ah,
    bool doStreaming=true, bool returnStreamedArray=false)
-> std::enable_if_t<decltype(has_advance(traj))::value || decltype(has_step(traj))::value,
                    trajectory::TemporalStreamedArray<typename std::decay_t<TRAJ>::StreamedArray>>
{
  if constexpr (decltype(has_step(traj))::value) { // it is of the Adaptive family
    if (p.dc) return run(std::forward<TRAJ>(traj),p.T,p.dc,p.sdf,p.ofn,p.initialFileName,p.precision,
                         p.streamInfo,p.firstStateStream,p.getParsedCommandLine(),
                         doStreaming,returnStreamedArray,std::forward<AutostopHandler>(ah));
    else if (!p.Dt) throw std::runtime_error("Nonzero dc or Dt required in trajectory::run");
  }
  if (!p.Dt) throw std::runtime_error("Nonzero Dt required in trajectory::run");
  if (p.NDt) return run(std::forward<TRAJ>(traj),p.NDt,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,
                        p.streamInfo,p.firstStateStream,p.getParsedCommandLine(),
                        doStreaming,returnStreamedArray,std::forward<AutostopHandler>(ah));
  else return run(std::forward<TRAJ>(traj),p.T,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,
                  p.streamInfo,p.firstStateStream,p.getParsedCommandLine(),
                  doStreaming,returnStreamedArray,std::forward<AutostopHandler>(ah));
}


template<typename TRAJ, typename ParsBase>
auto
run(TRAJ&& traj, const trajectory::Pars<ParsBase>& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return run(traj,p,trajectory::AutostopHandlerGeneric<typename std::decay_t<TRAJ>::StreamedArray>(p.autoStopEpsilon,p.autoStopEpsilonAbs,p.autoStopRepetition),
             doStreaming,returnStreamedArray);
}



namespace trajectory {


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


template <typename Trajectory>
struct ReadState
{
  static iarchive& _(Trajectory& traj, iarchive& iar)
  {
    SerializationMetadata sm;
    iar & sm;
    if (sm.trajectoryID==SerializationMetadata::ARRAY_ONLY) return traj.readFromArrayOnlyArchive(iar);
    else {
      if (sm != MakeSerializationMetadata<Trajectory>::_()) throw std::runtime_error("Trajectory mismatch in cppqedutils::trajectory::ReadState");
      return traj.stateIO(iar);
    }
  }
};


template <typename Trajectory>
struct WriteState
{
  static oarchive& _(Trajectory& traj, // cannot be const, because traj.stateIO is non-const
                     oarchive& oar)
  {
    return traj.stateIO(oar & MakeSerializationMetadata<Trajectory>::_());
  }  
};


template<typename Trajectory>
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


template<typename Trajectory>
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


/// Streams a limited set of relevant physical and numerical information about the actual state of Trajectory at the actual time instant
template <typename Trajectory>
auto stream(const Trajectory& traj, std::ostream& os,
            int precision ///< the precision (number of digits) of stream
            )
{
  const FormDouble fd(formdouble::positive(precision));
  auto res{traj.stream( os<<fd(traj.getTime())<<fd(traj.getDtDid()) , precision)};
  os<<std::endl; // Note: endl flushes the buffer
  return res;
}


/// Print header
template <typename Trajectory>
std::ostream& streamParameters(const Trajectory& traj, std::ostream& os)
{
  return traj.streamKey(traj.streamParameters(os)<<std::endl<<"Key to data:\nTrajectory\n 1. time\n 2. dtDid\n");
}


} } // cppqedutils::trajectory


#pragma GCC warning "TODO: This souldn’t call any TRAJ member function directly, only through traits classes"

template<typename TRAJ, typename LENGTH, typename DELTA, typename AutostopHandler>
auto
cppqedutils::trajectory::run(TRAJ&& traj, LENGTH length, DELTA streamFreq, unsigned stateStreamFreq,
                             const std::string& trajectoryFileName, const std::string& initialFileName,
                             int precision, bool streamInfo, bool firstStateStream,
                             const std::string& parsedCommandLine,
                             bool doStreaming, bool returnStreamedArray,
                             AutostopHandler&& autostopHandler) -> TemporalStreamedArray<typename std::decay_t<TRAJ>::StreamedArray>
{
  using namespace std;

  static constexpr bool
    endTimeMode=is_floating_point_v<LENGTH>,
    isDtMode=is_floating_point_v<DELTA>;
  
  TemporalStreamedArray<typename decay_t<TRAJ>::StreamedArray> res;
  
  ////////////////////////////////////////////////
  // Determining i/o streams, eventual state input
  ////////////////////////////////////////////////

  static const string stateExtension{".state"};
  const string stateFileName{trajectoryFileName+stateExtension};
  
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
    
  // this is formally a runtime expression, but since endTimeMode is constexpr, the compiler should be able to do this @ compile time
  const double timeToReach = (endTimeMode ? length : traj.getTime()+length*streamFreq);
  
  if (timeToReach && timeToReach<=traj.getTime()) return res;

  const shared_ptr<ostream> outstream{
    !streamToFile ?
    ( doStreaming ? 
      shared_ptr<ostream>(&cout,[](auto*){}) : // since cout is a system-wide object, this should be safe
      static_pointer_cast<ostream>(make_shared<boost::iostreams::stream<boost::iostreams::null_sink>>(boost::iostreams::null_sink{}) ) ) : 
    static_pointer_cast<ostream>(make_shared<ofstream>(trajectoryFileName.c_str(),ios_base::app))}; // regulates the deletion policy
  
  if (outstream->fail()) throw runtime_error("Trajectory stream opening error: "+trajectoryFileName);
  
  CommentingStream logStream{outstream};
  
  ostream& os=*outstream;
  IO_Manipulator::_(os);
  os<<setprecision(formdouble::actualPrecision(precision));

  CommentingStream commentingStream{outstream}; commentingStream<<setprecision(formdouble::actualPrecision(precision));
  
  ///////////////////////
  // Writing introduction
  ///////////////////////

  if (streamInfo) {
    if (!continuing) {
      if (parsedCommandLine!="") commentingStream<<parsedCommandLine<<endl<<endl;
      streamParameters(traj,commentingStream<<versionHelper())
        <<endl<<"Run Trajectory up to time "<<timeToReach
        <<" -- Stream period: "<<streamFreq<< (isDtMode ? "" : " timestep") <<endl<<endl;
    }
    else
      commentingStream<<"Continuing from time "<<traj.getTime()<<" up to time "<<timeToReach<<endl;
  }

  if (!timeToReach) {stream(traj,os,precision); return res;}

  //////////////////////////////
  // Mid section: the actual run
  //////////////////////////////

  const std::shared_ptr<ostream> ofs = !streamToFile ? std::make_shared<ofstream>() : openStateFileWriting(stateFileName);

  bool
    stateSaved=false,  // signifies whether the state has already been saved for the actual time instant of the trajectory
    arrayStreamed=false; // signifies whether the expectation values have already been streamed ”

  try {

    for (long count=0, stateCount=0; (endTimeMode ? traj.getTime()<length : count<=length); ++count) {

      if (count) {
        // advance trajectory
        if constexpr (!isDtMode) {traj.step(length-traj.getTime(),logStream);}
        else {
          if constexpr (endTimeMode) advance(traj,std::min(streamFreq,length-traj.getTime()),logStream);
          else advance(traj,streamFreq,logStream);
        }
        stateSaved=arrayStreamed=false;
      }

      if (!count || [&]() {
        if constexpr (isDtMode) return true;
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
          auto streamReturn{stream(traj,os,precision)};
          if (returnStreamedArray) res.emplace_back(traj.getTime(),traj.getDtDid(),get<1>(streamReturn));
          autostopHandler(get<1>(streamReturn));
        }
      }
    } // end of main for loop

  } catch (const StoppingCriterionReachedException& except) {commentingStream<<"Stopping criterion has been reached"<<endl;}
  
  if (!arrayStreamed) stream(traj,os,precision); // Stream at the end instant if stream has not happened yet

  //////////////////////////////////////////
  // Logging on end, saving trajectory state
  //////////////////////////////////////////
  
  traj.logOnEnd(commentingStream);
  if (!stateSaved) writeViaSStream(traj,ofs);
  
  return res;
  
}



#endif // CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED
