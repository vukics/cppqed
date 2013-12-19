/// \briefFile{Defines the basic classes of the trajectory-bundle}
// -*- C++ -*-
#ifndef UTILS_TRAJECTORY_H_INCLUDED
#define UTILS_TRAJECTORY_H_INCLUDED

#include "TrajectoryFwd.h"

#include "Archive.h"
#include "Exception.h"
#include "Evolved.h"

#include <iosfwd>
#include <string>


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
         bool firstStateDisplay ///< governs whether the state is displayed at time zero (important if \link Trajectory::writeState state display\endlink is costly)
        );

void run(Trajectory &, long   nDt , double deltaT, unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay);

template<typename A>
void run(Adaptive<A>&, double time, int dc       , unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay);

void run(Trajectory &, const ParsRun&);

template<typename A>
void run(Adaptive<A>&, const ParsRun&);



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


template<typename T>
void writeViaSStream(const T&, std::ofstream*);
template<typename T>
void  readViaSStream(      T&, std::ifstream&);
SerializationMetadata readMeta(std::ifstream&);


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


inline double initialTimeStep(double highestFrequency) {return 1./(10.*highestFrequency);}
// A heuristic determination of the inital timestep from the highest frequency of a physical system.



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
  
};


///////////
//
// Adaptive
//
///////////

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
   * The dimensions check can be circumvented by setting the trajectoryID to
   * SerializationMetadata::ArrayOnly in the EvolvedIO's metadata. This is done for example
   * in the python I/O interface, because when reading in a state in python we
   * generally have no idea about the dimensions.
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


/// Adaptive is more than Evolved only in that it takes into account the need for communicating towards the user from time to time during propagation
template<typename A>
class Adaptive : public trajectory::AdaptiveIO<A>, public virtual Trajectory
{
public:
  // needed to break ambiguity with identically named functions in AdaptiveIO
  using Trajectory::readState; using Trajectory::writeState; using Trajectory::getTime;

  typedef evolved::Evolved<A> Evolved;

  void step(double deltaT) {step_v(deltaT);}
  
  virtual ~Adaptive() {}

protected:
  using AdaptiveIO<A>::meta_;

  Adaptive(A&, typename Evolved::Derivs, double, double, double    , const A&, const evolved::Maker<A>&);

  Adaptive(A&, typename Evolved::Derivs, double, const ParsEvolved&, const A&, const evolved::Maker<A>&);

  typedef typename Evolved::ConstPtr ConstPtr;
  typedef typename Evolved::     Ptr      Ptr;
  
  const ConstPtr getEvolved() const {return ConstPtr(evolved_);}
  const      Ptr getEvolved()       {return          evolved_ ;}

  double getDtTry() const {return evolved_->getDtTry();}

  std::ostream& displayParameters_v(std::ostream&) const;

  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       final;
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const final;

  const std::string trajectoryID() const  {return trajectoryID_v();}

  virtual cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar)       {return iar;}
  virtual cpputils::oarchive& writeStateMore_v(cpputils::oarchive &oar) const {return oar;}

private:

  double getDtDid_v() const {return evolved_->getDtDid();}

  void evolve_v(double deltaT) {evolved::evolve<Adaptive>(*this,deltaT);}

  double getTime_v() const {return evolved_->getTime();}

  virtual void step_v(double deltaT) = 0;

  virtual const std::string trajectoryID_v() const = 0;

  const Ptr evolved_;

  const double dtInit_;

};


} // trajectory


#endif // UTILS_TRAJECTORY_H_INCLUDED
