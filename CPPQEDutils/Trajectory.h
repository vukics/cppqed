// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the basic classes of the trajectory-bundle}
#ifndef CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED

#include "Archive.h"
#include "CommentingStream.h"
#include "Evolved.h"
#include "FormDouble.h"
#include "ParsTrajectory.h"

#include <boost/hana.hpp>

#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>

namespace hana=boost::hana;


/// The trajectory-bundle
namespace trajectory {

/// Type-erased aggregate of information about a trajectory-state archive \see AdaptiveIO
/**
 * \par Rationale
 * Each archive is self-contained, so that it contains its own metadata.
 * Then, the state files can eventually even be split into several archives, each containing a single self-contained state.
 */
struct SerializationMetadata
{
  SerializationMetadata(std::string type=UNSPECIFIED, std::string id=UNSPECIFIED, int r=0)
    : protocolVersion(0),
      rank(r),
      typeID(type),
      trajectoryID(id)
  {};

  int protocolVersion;
  int rank;
  std::string typeID;
  std::string trajectoryID;

  static const std::string UNSPECIFIED;
  static const std::string ARRAY_ONLY;

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & protocolVersion & rank & typeID & trajectoryID;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

};

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

template<typename T>
void writeViaSStream(const T&, std::shared_ptr<std::ostream>);
template<typename T>
void  readViaSStream(      T&, std::shared_ptr<std::istream>);
SerializationMetadata readMeta(std::shared_ptr<std::istream>); ///< Needed separately for the Python i/o


struct StoppingCriterionReachedException {};


/// A heuristic determination of the inital timestep from the highest frequency of a physical system.
inline double initialTimeStep(double highestFrequency) {return 1./(10.*highestFrequency);}


/// The base class of the trajectory-bundle condensing the quintessential characteristics of any trajectory
/**
 * \see run(Trajectory&, double, double, unsigned, const std::string&, const std::string&, int, bool, bool) which manifests all the principal functionalities of this class
 *
 * The i/o of the entire state is performed via \refBoost{Boost.Serialization,serialization}, and is disabled if this library is not installed
 *
 * \todo Consider the possibility of cloning Trajectories
 */
template<typename SA>
class Trajectory
{
protected:
  Trajectory() = default;
  Trajectory(const Trajectory&) = delete; Trajectory& operator=(const Trajectory&) = delete;
  Trajectory(Trajectory&&) = default; Trajectory& operator=(Trajectory&&) = default;
  
public:
  /// Propagation for a time interval of exactly deltaT
  /** \note logStream must always be flushed (e.g. with `std::endl`) before anything else gets to write on the same underlying stream  */
  void evolve(double deltaT, std::ostream& logStream) {evolve_v(deltaT, logStream);}

  using StreamedArray=SA;
  using StreamReturnType=std::tuple<std::ostream&,SA>;
  
  /// Streams a limited set of relevant physical and numerical information about the actual state of Trajectory at the actual time instant
  StreamReturnType stream(std::ostream& os, int precision ///< the precision (number of digits) of stream
                         ) const
  {
    const FormDouble fd(formdouble::positive(precision));
    auto res{stream_v( os<<fd(getTime())<<fd(getDtDid()) , precision)};
    std::get<0>(res)<<std::endl;
    return res; // Note: endl flushes the buffer
  }


  /// \name Getters
  //@{
  double getTime() const {return getTime_v();} ///< actual time instant
  double getDtDid() const {return getDtDid_v();} ///< last perfomed timestep
  //@}

  std::ostream& streamParameters(std::ostream& os) const ///< print header
  {
    size_t i=3;
    return streamKey_v(streamParameters_v(os)<<std::endl<<"Key to data:\nTrajectory\n 1. time\n 2. dtDid\n" , i);
  }
  
  std::ostream& logOnEnd(std::ostream& os) const {return logOnEnd_v(os);} ///< print a log at the end summarizing overall (e.g. time-averaged) physical and numerical data during the run

  /// \name Entire state i/o
  //@{
  cpputils::iarchive&  readState(cpputils::iarchive& iar)       {return  readState_v(iar);} ///< read from an archive
  cpputils::oarchive& writeState(cpputils::oarchive& oar) const {return writeState_v(oar);} ///< write to an archive
  //@}

  virtual ~Trajectory() {}

private:
  virtual void evolve_v(double, std::ostream&) = 0; 
  virtual double getTime_v() const = 0;
  virtual double getDtDid_v() const = 0;
  
  virtual std::ostream& streamParameters_v(std::ostream&) const = 0;
  virtual StreamReturnType stream_v(std::ostream&, int) const = 0;
  virtual std::ostream& streamKey_v(std::ostream&, size_t&) const = 0;

  virtual cpputils::iarchive&  readState_v(cpputils::iarchive&) = 0;
  virtual cpputils::oarchive& writeState_v(cpputils::oarchive&) const = 0;

  virtual std::ostream& logOnEnd_v(std::ostream& os) const {return os;}

};


///////////
//
// Adaptive
//
///////////

/// Corresponds to evolved::EvolvedIO, adding the capability of state serialization involving a SerializationMetadata instant
/** This class exists for the sake of Python I/O */
template<typename A>
class AdaptiveIO
{
protected:
  AdaptiveIO(AdaptiveIO&&) = default; AdaptiveIO& operator=(AdaptiveIO&&) = default;
  
private:
  typedef evolved::EvolvedIO<A> EvolvedIO;

  typedef typename EvolvedIO::Ptr Ptr;

public:
  AdaptiveIO(Ptr);

  /// Read in the EvolvedIO from a cpputils::iarchive.
  /**
   * Two conformity checks are performed for to the array we try to read in:
   * - whether the rank matches (throws RankMismatchException if not)
   * - whether the dimensions match (throws DimensionsMismatchException if not)
   * The dimensions check can be circumvented by setting the trajectoryID to SerializationMetadata::ArrayOnly in the EvolvedIO's metadata.
   * This is done for example in the python I/O interface, because when reading in a state in python we generally have no idea about the dimensions.
   */
  cpputils::iarchive&  readState(cpputils::iarchive& iar);
  /// Write the EvolvedIO to a cpputils::oarchive
  cpputils::oarchive& writeState(cpputils::oarchive& oar) const;

  /// Returns the time of the underlying EvolvedIO
  double getTime() const {return evolvedIO_->getTime();}

protected:
  const Ptr getEvolvedIO() const {return evolvedIO_;} ///< note: not the same const-correctness as in Adaptive

  mutable SerializationMetadata meta_;

private:
  const Ptr evolvedIO_;

};



/// Adaptive is basically an evolved::Evolved wrapped into the Trajectory interface
/** The class stores an evolved::Evolved instance by shared pointer
 * 
 * \tparam BASE either Trajectory or Averageable (cf. StochasticTrajectory.h)
 * 
 */
template<typename A, typename BASE>
class Adaptive : public trajectory::AdaptiveIO<A>, public BASE
{
protected:
  Adaptive(Adaptive&&) = default; Adaptive& operator=(Adaptive&&) = default;

public:
  // needed to break ambiguity with identically named functions in AdaptiveIO
  using BASE::readState; using BASE::writeState; using BASE::getTime;

  typedef evolved::Evolved<A> Evolved;

  /// corresponding to Evolved::step, it takes a single adaptive step
  /** It does not delegate directly to Evolved::step, as usually trajectories need to do more for a step than just propagating the ODE: instead, it is kept purely virtual */
  void step(double deltaT, ///< *maximum* length of the timestep
            std::ostream& logStream);
  
  virtual ~Adaptive() {}

protected:
  using AdaptiveIO<A>::meta_;

  typedef typename Evolved::Derivs Derivs;
  
  /// Constructor taking the same parameters as needed to operate evolved::Maker
  template<typename ARRAY, typename... BaseInitializationPack>
  Adaptive(ARRAY&& y, Derivs derivs, double dtInit, int logLevel, double epsRel, double epsAbs,
           const A& scaleAbs, const evolved::Maker<A>& maker,
           BaseInitializationPack&&... bip)
    : AdaptiveIO<A>{maker(std::forward<ARRAY>(y),derivs,dtInit,epsRel,epsAbs,scaleAbs)},
      BASE{std::forward<BaseInitializationPack>(bip)...},
      evolved_{std::dynamic_pointer_cast<Evolved>(AdaptiveIO<A>::getEvolvedIO())},
      dtInit_(dtInit), logLevel_(logLevel) {}

  typedef typename Evolved::ConstPtr ConstPtr;
  typedef typename Evolved::     Ptr      Ptr;

  /// \name Getters
  //@{
  const ConstPtr getEvolved() const {return ConstPtr(evolved_);}
  const      Ptr getEvolved()       {return          evolved_ ;}

  double getDtTry() const {return evolved_->getDtTry();}
  //@}
  
  /// redirected to a pure virtual, this is needed for \link SerializationMetadata serialization of trajectory metadata\endlink
  const std::string trajectoryID() const  {return trajectoryID_v();}

  std::ostream& streamParameters_v(std::ostream&) const override;

  virtual cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar)       {return iar;} ///< hook into Trajectory::readState
  virtual cpputils::oarchive& writeStateMore_v(cpputils::oarchive &oar) const {return oar;} ///< hook into Trajectory::writeState

private:
  /// calls AdaptiveIO::readState, checks post-conditions in metadata, calls readStateMore_v, and finally sets `dtTry` in `evolved_`
  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       final;
  /// calls AdaptiveIO::writeState followed by writeStateMore_v
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const final;

  double getDtDid_v() const final {return evolved_->getDtDid();}

  void evolve_v(double deltaT, std::ostream& logStream) final {evolved::evolve<Adaptive>(*this,deltaT,logStream);}

  double getTime_v() const final {return evolved_->getTime();}

  virtual void step_v(double deltaT, std::ostream& logStream) = 0;

  virtual const std::string trajectoryID_v() const = 0;
  
  std::ostream& logOnEnd_v(std::ostream& os) const final {if (logLevel_) evolved_->logOnEnd(os); return logMoreOnEnd_v(os);}

  virtual std::ostream& logMoreOnEnd_v(std::ostream& os) const {return os;} ///< hook into Trajectory::logOnEnd
  
  const Ptr evolved_;

  const double dtInit_;
  
  const int logLevel_;

};


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
  AutostopHandlerGeneric(double autoStopEpsilon, unsigned autoStopRepetition) : autoStopEpsilon_{autoStopEpsilon}, autoStopRepetition_{autoStopRepetition}, averages_{}, queue_{} {}

  void operator()(const SA& streamedArray)
  {
    if (!autoStopRepetition_) return;
    
    if (!averages_.size()) { // This means that no stream has yet occured: the size of averages_ must be determined
      averages_.resize(streamedArray.size());
      averages_=0.;
    }
    else {
      if (queue_.size()==autoStopRepetition_) {
        if (max(abs(averages_-streamedArray)/(abs(averages_)+abs(streamedArray)))<autoStopEpsilon_) throw StoppingCriterionReachedException();
        averages_=averages_+(streamedArray-queue_.front())/double(autoStopRepetition_); // update the averages set for next step
        queue_.pop(); // remove obsolate first element of the queue
      }
      else averages_=(double(queue_.size())*averages_+streamedArray)/double(queue_.size()+1); // build the initial averages set to compare against

      queue_.push(streamedArray); // place the new item to the back of the queue

    }

  }

private:
  const double autoStopEpsilon_;
  const unsigned autoStopRepetition_;

  SA averages_;
  std::queue<SA> queue_;

};


template<typename SA>
struct AutostopHandlerNoOp
{
  void operator()(const SA&) {}
};


template<typename SA>
using TemporalStreamedArray=std::list<std::tuple<double,double,SA>>;


/// Running in deltaT mode (streams in equal time intervals) for a certain time \related Trajectory
/**
 * This function manifests all the basic features of Adaptive and the whole idea behind the trajectory bundle.
 *
 * A Trajectory can
 * - be \link Trajectory::evolve evolved\endlink (propagated in time by given time intervals)
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
template<typename T, // type of trajectory
         typename L, // type specifying the length of the run
         typename D, // type specifying the frequency of stream
         typename AutostopHandler // should support operator()(const typename T::TemporalStreamedArray &)
         >
TemporalStreamedArray<typename std::decay_t<T>::StreamedArray>
run(T&& traj, ///< the trajectory to run
    L length, ///< length of run
    D streamFreq, ///< interval between two \link Trajectory::stream streamings\endlink
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


inline auto has_evolve = hana::is_valid([](auto&& obj) -> decltype(obj.evolve(1.0,std::clog)) { });
inline auto has_step = hana::is_valid([](auto&& obj) -> decltype(obj.step(1.0,std::clog)) { });

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
template<typename TRAJ>
auto
run(TRAJ&& traj, const ParsRun& p,
    bool doStreaming=true, bool returnStreamedArray=false)
-> std::enable_if_t<decltype(has_evolve(traj))::value,TemporalStreamedArray<typename std::decay_t<TRAJ>::StreamedArray>>
{
  using AHG=AutostopHandlerGeneric<typename std::decay_t<TRAJ>::StreamedArray>;
  
  if constexpr (decltype(has_step(traj))::value) { // it is of the Adaptive family
    if (p.dc) return run(std::forward<TRAJ>(traj),p.T,p.dc,p.sdf,p.ofn,p.initialFileName,p.precision,
                         p.streamInfo,p.firstStateStream,p.getParsedCommandLine(),
                         doStreaming,returnStreamedArray,
                         AHG(p.autoStopEpsilon,p.autoStopRepetition));
    else if (!p.Dt) throw std::runtime_error("Nonzero dc or Dt required in trajectory::run");
  }
  if (!p.Dt) throw std::runtime_error("Nonzero Dt required in trajectory::run");
  if (p.NDt) return run(std::forward<TRAJ>(traj),p.NDt,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,
                        p.streamInfo,p.firstStateStream,p.getParsedCommandLine(),
                        doStreaming,returnStreamedArray,
                        AHG(p.autoStopEpsilon,p.autoStopRepetition));
  else return run(std::forward<TRAJ>(traj),p.T,p.Dt,p.sdf,p.ofn,p.initialFileName,p.precision,
                  p.streamInfo,p.firstStateStream,p.getParsedCommandLine(),
                  doStreaming,returnStreamedArray,
                  AHG(p.autoStopEpsilon,p.autoStopRepetition));
}


} // trajectory



#endif // CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED
