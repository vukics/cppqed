// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the basic classes of the trajectory-bundle}
#ifndef CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED

#include "TrajectoryFwd.h"

#include "Archive.h"
#include "CommentingStream.h"
#include "Exception.h"
#include "Evolved.h"
#include "SmartPtr.h"

#include <iosfwd>
#include <string>
#include <boost/shared_ptr.hpp>


/// The trajectory-bundle
namespace trajectory {


/// Running in deltaT mode (displays in equal time intervals) for a certain time \related Trajectory
/**
 * This function manifests all the basic features of Adaptive and the whole idea behind the trajectory bundle.
 *
 * A Trajectory can
 * - be \link Trajectory::evolve evolved\endlink (propagated in time by given time intervals)
 * - \link Trajectory::readState perform i/o of its entire state\endlink, that is, a bunch of information necessary for resuming a Trajectory from a certain time instant
 * - \link Trajectory::display display relevant physical and numerical information\endlink about its actual state at any time (e.g. a set of quantum averages in the case of a quantum trajectory)
 *
 * \note While the entire state can be huge (e.g. the state vector or density operator in the case of a quantum trajectory) the relevant information in an actual numerical experiment
 * is usually much less (a set of quantum averages that doesn’t entirely define the state).
 *
 * Furthermore, a Trajectory can
 * - provide information about its \link Trajectory::getTime time\endlink and \link Trajectory::getDtDid last performed timestep\endlink
 * - \link Trajectory::displayParameters print a header\endlink summarizing its physical and numerical parameters together with a key to the set of relevant physical information displayed during the run
 * - \link Trajectory::logOnEnd print a log\endlink at the end summarizing overall (e.g. time-averaged) physical and numerical data during the run
 *
 * \see Simulated for a full generic implementation of Trajectory together with a small tutorial
 * 
 * \todo Consider taking Trajectory by rvalue reference
 *
 */
void run(Trajectory & trajectory, ///< the trajectory to run
         double time, ///< end time
         double deltaT, ///< time interval between two \link Trajectory::display displays\endlink
         unsigned sdf, ///< number of \link Trajectory::display displays\endlink between two \link Trajectory::writeState state displays\endlink
         const std::string& ofn, ///< name of the output file for \link Trajectory::display displays\endlink — if empty, display to standard output; \link Trajectory::writeState state displays\endlink into file named `ofn.state`
         const std::string& initialFileName, ///< name of file containing initial condition state for the run
         int precision, ///< governs the overall precision (number of digits) of outputs in \link Trajectory::display displays\endlink
         bool displayInfo, ///< governs whether a \link Trajectory::displayParameters header\endlink is displayed at the top of the output
         bool firstStateDisplay, ///< governs whether the state is displayed at time zero (important if \link Trajectory::writeState state display\endlink is costly)
         double autoStopEpsilon, ///< relative precision for autostopping
         unsigned autoStopRepetition, ///< number of displayed lines repeated within relative precision before autostopping – 0 means no autostopping
         const std::string& parsedCommandLine
        );


/// Same as \link Trajectory::run above\endlink but runs for a certain number of time intervals deltaT \related Trajectory
/**
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
 */
void run(Trajectory &, long nDt, ///< the end time of the trajectory will be nDt*deltaT
         double deltaT, unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay,
         double autoStopEpsilon, unsigned autoStopRepetition, const std::string& parsedCommandLine);


/// Another version of \link Trajectory::run `run`\endlink for running in dc-mode \related Adaptive
/**
 * Since in addition to the Trajectory interface, Adaptive has also the capability to be propagated over a \link Adaptive::step single adaptive timestep\endlink, it is possible to count the
 * individual ODE steps. Running in dc-mode means that a fixed number of adaptive steps are performed between each display. Hence, we get denser displays in time when the timestep is small,
 * that is, when the important things are happening in the dynamics.
 * 
 */
template<typename A, typename BASE>
void run(Adaptive<A,BASE>&, double time, int dc, ///< number of adaptive timesteps taken between two displays
         unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay,
         double autoStopEpsilon, unsigned autoStopRepetition, const std::string& parsedCommandLine);


/// Dispatcher \related Trajectory
/**
 * Runs
 * - run(Trajectory&, long, double, unsigned, const std::string&, const std::string&, int, bool, bool) if `p.NDt` is nonzero and
 * - run(Trajectory&, double, double, unsigned, const std::string&, const std::string&, int, bool, bool) otherwise
 * 
 * \note This means that ParsRun::NDt takes precedence over ParsRun::T
 *
 */
void run(Trajectory &, const ParsRun& p);


/// Dispatcher \related Adaptive
/**
 * - Runs run(Adaptive&, double, int, unsigned, const std::string&, const std::string&, int, bool, bool) if `p.dc` is nonzero (dc-mode)
 * - delegates to run(Trajectory&, const ParsRun& p) otherwise (deltaT-mode)
 * 
 * \note This means that ParsRun::dc takes precedence over ParsRun::Dt
 * 
 */
template<typename A, typename BASE>
void run(Adaptive<A,BASE>&, const ParsRun&);


/// Aggregate of information about a trajectory-state archive \see AdaptiveIO
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
#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & protocolVersion & rank & typeID & trajectoryID;}
#endif // DO_NOT_USE_BOOST_SERIALIZATION

};

/// Opens a file with the given filename for reading and returns an associated std::istream object.
/**
 * If C++QED is compiled with compression support (boost iostreams library needed), then
 * the statefile can be in bzip2-format, the returned istream will handle this automatically.
 *
 * If the file cannot be opened, StateFileOpeningException is raised.
 */
boost::shared_ptr<std::istream> openStateFileReading(const std::string &filename);

/// Opens a file with the given filename for writing (appending) and returns an associated std::ostream object.
/**
 * If C++QED is compiled with compression support (boost iostreams library needed), then
 * the statefile can be in bzip2-format, the returned ostream will handle this automatically.
 * If the file does not yet exist, bzip2 is used if supported.
 *
 * If the file cannot be opened, StateFileOpeningException is raised.
 */
boost::shared_ptr<std::ostream> openStateFileWriting(const std::string &filename, const std::ios_base::openmode mode=std::ios_base::app | std::ios_base::binary);

template<typename T>
void writeViaSStream(const T&, std::ostream*);
template<typename T>
void  readViaSStream(      T&, std::istream*);
SerializationMetadata readMeta(std::istream*); ///< Needed separately for the Python i/o


class StoppingCriterionReachedException : public cpputils::Exception {};

/// Raised when the rank of a trajectory we try to read in from file does not match
class RankMismatchException : public cpputils::Exception {};

/// Raised when the dimensions of a trajectory state we try to read in from file does not match
class DimensionsMismatchException : public cpputils::Exception {};

/// Raised when the trajectory type we try to read in from file does not match.
class TrajectoryMismatchException : public cpputils::Exception {};

class TrajectoryFileOpeningException : public cpputils::TaggedException
{
public:
  TrajectoryFileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

};


class StateFileOpeningException : public cpputils::TaggedException
{
public:
  StateFileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

};


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
class Trajectory : private boost::noncopyable
{
public:
  class CommentingStreamUnsetException : public cpputils::Exception {};
  
  /// Propagation for a time interval of exactly deltaT
  void evolve(double deltaT) {evolve_v(deltaT);}

  /// Displays a limited set of relevant physical and numerical information about the actual state of Trajectory at the actual time instant
  std::ostream& display(std::ostream&, int precision ///< the precision (number of digits) of display
                       ) const;

  /// \name Getters
  //@{
  double getTime() const {return getTime_v();} ///< actual time instant
  double getDtDid() const {return getDtDid_v();} ///< last perfomed timestep
  //@}

  std::ostream& displayParameters(std::ostream& os) const; ///< print header

  std::ostream& logOnEnd(std::ostream& os) const {return logOnEnd_v(os);} ///< print a log at the end summarizing overall (e.g. time-averaged) physical and numerical data during the run

  /// \name Entire state i/o
  //@{
  cpputils::iarchive&  readState(cpputils::iarchive& iar)       {return  readState_v(iar);} ///< read from an archive
  cpputils::oarchive& writeState(cpputils::oarchive& oar) const {return writeState_v(oar);} ///< write to an archive
  //@}

  virtual ~Trajectory() {}

  void setLogStreamDuringRun(boost::shared_ptr<std::ostream> os) {commentingStream_.p=os; commentingStream_.s=*os;}

protected:
  /// The stream hence obtained must always be flushed (e.g. with `std::endl`) before anything else gets to write on the same underlying stream
  std::ostream& getLogStreamDuringRun() const {return commentingStream_.s;}
  
private:
  virtual void            evolve_v(double)       = 0; // A step of exactly deltaT
  virtual double         getTime_v()       const = 0;
  virtual double        getDtDid_v()       const = 0;
  
  virtual std::ostream& displayParameters_v(std::ostream&         ) const = 0;
  virtual std::ostream& display_v          (std::ostream&, int    ) const = 0;
  virtual std::ostream& displayKey_v       (std::ostream&, size_t&) const = 0;

  virtual cpputils::iarchive&  readState_v(cpputils::iarchive&) = 0;
  virtual cpputils::oarchive& writeState_v(cpputils::oarchive&) const = 0;

  virtual std::ostream& logOnEnd_v(std::ostream& os) const {return os;}

  struct {
    boost::shared_ptr<std::ostream> p;
    mutable cpputils::CommentingStream s{std::clog};
  } commentingStream_;

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
   * The dimensions check can be circumvented by setting the trajectoryID to SerializationMetadata::ArrayOnly in the EvolvedIO's metadata. This is done for example in the python I/O interface, 
   * because when reading in a state in python we generally have no idea about the dimensions.
   * 
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
public:
  // needed to break ambiguity with identically named functions in AdaptiveIO
  using BASE::readState; using BASE::writeState; using BASE::getTime;

  typedef evolved::Evolved<A> Evolved;

  /// corresponding to Evolved::step, it takes a single adaptive step
  /** It does not delegate directly to Evolved::step, as usually trajectories need to do more for a step than just propagating the ODE: instead, it is kept purely virtual */
  void step(double deltaT ///< *maximum* length of the timestep
           );
  
  virtual ~Adaptive() {}

protected:
  using AdaptiveIO<A>::meta_;

  typedef typename Evolved::Derivs Derivs;
  
  /// Constructor taking the same parameters as needed to operate evolved::Maker
  template<typename... BaseInitializationPack>
  Adaptive(A&, Derivs, double, int, double, double, const A&, const evolved::Maker<A>&, BaseInitializationPack&&...);

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

  std::ostream& displayParameters_v(std::ostream&) const override;

  virtual cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar)       {return iar;} ///< hook into Trajectory::readState
  virtual cpputils::oarchive& writeStateMore_v(cpputils::oarchive &oar) const {return oar;} ///< hook into Trajectory::writeState

private:
  /// calls AdaptiveIO::readState, checks post-conditions in metadata, calls readStateMore_v, and finally sets `dtTry` in `evolved_`
  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       final;
  /// calls AdaptiveIO::writeState followed by writeStateMore_v
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const final;

  double getDtDid_v() const final {return evolved_->getDtDid();}

  void evolve_v(double deltaT) final {evolved::evolve<Adaptive>(*this,deltaT);}

  double getTime_v() const final {return evolved_->getTime();}

  virtual void step_v(double deltaT) = 0;

  virtual const std::string trajectoryID_v() const = 0;
  
  std::ostream& logOnEnd_v(std::ostream& os) const final {if (logLevel_) evolved_->logOnEnd(os); return logMoreOnEnd_v(os);}

  virtual std::ostream& logMoreOnEnd_v(std::ostream& os) const {return os;} ///< hook into Trajectory::logOnEnd
  
  const Ptr evolved_;

  const double dtInit_;
  
  const int logLevel_;

};


} // trajectory


#endif // CPPQEDCORE_UTILS_TRAJECTORY_H_INCLUDED
