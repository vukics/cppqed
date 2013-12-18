// -*- C++ -*-
// Adaptive is more than Evolved only in that it takes into account the need for communicating towards the user from time to time during Evolution.
// This is manifested in the abstract function Display.
// Consider copying (cloning) of Trajectories

#ifndef UTILS_TRAJECTORY_H_INCLUDED
#define UTILS_TRAJECTORY_H_INCLUDED

#include "TrajectoryFwd.h"

#include "Archive.h"
#include "Exception.h"
#include "Evolved.h"

#include <iosfwd>
#include <string>



namespace trajectory {


void run(Trajectory &, double time, double deltaT, unsigned sdf, const std::string& ofn, const std::string& initialFileName, int precision, bool displayInfo, bool firstStateDisplay);

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


/////////////
//
// Trajectory
//
/////////////


class Trajectory : private boost::noncopyable
{
public:
  std::ostream& display   (std::ostream&, int) const;

  void evolve(double deltaT) {evolve_v(deltaT);} // A step of exactly deltaT

  double getTime() const {return getTime_v();}

  double getDtDid() const {return getDtDid_v();}
  std::ostream& displayParameters(std::ostream& os) const;

  
  std::ostream& logOnEnd(std::ostream& os) const {return logOnEnd_v(os);}
  
  cpputils::iarchive&  readState(cpputils::iarchive& iar)       {return  readState_v(iar);}
  cpputils::oarchive& writeState(cpputils::oarchive& oar) const {return writeState_v(oar);};
  
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
  void resetInitialDtTry() {evolved_->setDtTry(dtInit_);}

  std::ostream& displayParameters_v(std::ostream&) const;

  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       final;
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const final;

  const std::string trajectoryID() const  {return trajectoryID_v();}

  virtual cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar)       {return iar;}
  virtual cpputils::oarchive& writeStateMore_v(cpputils::oarchive &oar) const {return oar;}

private:
  using AdaptiveIO<A>::getEvolvedIO;

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
