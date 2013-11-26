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

#include <boost/utility.hpp>

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
  SerializationMetadata(std::string id=UNSPECIFIED)
    : protocolVersion(0),
      rank(0),
      trajectoryID(id)
  {};
  int protocolVersion;
  int rank;
  std::string trajectoryID;
  
  static const std::string UNSPECIFIED;
  static const std::string ARRAY_ONLY;

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & protocolVersion & rank & trajectoryID;}
#endif // DO_NOT_USE_BOOST_SERIALIZATION
};

namespace details
{

void writeNextArchive(std::ofstream*, const std::ostringstream&);
void readNextArchive(std::ifstream&, std::istringstream&);

} // details

void writeViaSStream(const Trajectory&, std::ofstream*);
void  readViaSStream(      Trajectory&, std::ifstream&);
SerializationMetadata readMeta(std::ifstream&);


class StoppingCriterionReachedException : public cpputils::Exception {};


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
  std::ostream& displayKey(std::ostream&     ) const;

  void evolve(double deltaT) {evolve_v(deltaT);} // A step of exactly deltaT

  double getTime() const {return getTime_v();}

  double getDtDid() const {return getDtDid_v();}

  std::ostream& displayParameters(std::ostream& os) const {return displayParameters_v(os);}
  
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
public:
  typedef evolved::EvolvedIO<A>   EvolvedIO;
  AdaptiveIO(typename EvolvedIO::Ptr evolvedIO) 
    : meta_(SerializationMetadata::ARRAY_ONLY), evolvedIO_(evolvedIO) {};

  typedef typename EvolvedIO::ConstPtr ConstPtr;
  typedef typename EvolvedIO::Ptr Ptr;

  const ConstPtr getEvolvedIO() const {return ConstPtr(evolvedIO_);}
  const      Ptr getEvolvedIO()       {return          evolvedIO_ ;}

  cpputils::iarchive&  readArrayState(cpputils::iarchive& iar)       {return iar & meta_ & *evolvedIO_;}
  cpputils::oarchive& writeArrayState(cpputils::oarchive& oar) const {return oar & meta_ & *evolvedIO_;}
protected:
  mutable SerializationMetadata meta_;
private:
  Ptr evolvedIO_;
};

template<typename A>
using EvolvedPtrBASE = boost::base_from_member<typename evolved::Evolved<A>::Ptr>;

template<typename A>
class Adaptive : private EvolvedPtrBASE<A>, public trajectory::AdaptiveIO<A>, public virtual Trajectory
{
public:
  // Some parameter-independent code could still be factored out, but probably very little
  
  typedef evolved::Evolved<A>       Evolved;
  // typedef trajectory::AdaptiveIO<A> AdaptiveIO;

  void step(double deltaT) {step_v(deltaT);}
  
  virtual ~Adaptive() {}

protected:
  using AdaptiveIO<A>::meta_;

  Adaptive(A&, typename Evolved::Derivs, double, double, double    , const A&, const evolved::Maker<A>&);

  Adaptive(A&, typename Evolved::Derivs, double, const ParsEvolved&, const A&, const evolved::Maker<A>&);

  typedef typename Evolved::ConstPtr ConstPtr;
  typedef typename Evolved::Ptr Ptr;
  
  const ConstPtr getEvolved() const {return ConstPtr(evolved_);}
  const      Ptr getEvolved()       {return          evolved_ ;}

  double getDtTry() const {return evolved_->getDtTry();}

  std::ostream& displayParameters_v(std::ostream&) const;

  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       final;
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const final;

  std::string trajectoryID() const  {return trajectoryID_v();}

private:
  using AdaptiveIO<A>::readArrayState;
  using AdaptiveIO<A>::writeArrayState;

  typedef EvolvedPtrBASE<A> EvolvedPtrBase;

  double getDtDid_v() const {return evolved_->getDtDid();}

  void evolve_v(double deltaT) {evolved::evolve<Adaptive>(*this,deltaT);}

  double getTime_v() const {return evolved_->getTime();}

  virtual void step_v(double deltaT) = 0;
  // Prefer purely virtual functions, so that there is no danger of forgetting to override them. Very few examples anyway for a trajectory wanting to perform only a step of Evolved.
  // (Only Simulated, but neither Master, nor MCWF_Trajectory)

  virtual cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar)       {return iar;}
  virtual cpputils::oarchive& writeStateMore_v(cpputils::oarchive &oar) const {return oar;}

  virtual std::string trajectoryID_v() const = 0;

  const Ptr evolved_;

};


} // trajectory


#endif // UTILS_TRAJECTORY_H_INCLUDED
