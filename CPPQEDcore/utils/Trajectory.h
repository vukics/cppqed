// -*- C++ -*-
// Adaptive is more than Evolved only in that it takes into account the need for communicating towards the user from time to time during Evolution.
// This is manifested in the abstract function Display.
// Consider copying (cloning) of Trajectories

#ifndef UTILS_INCLUDE_TRAJECTORY_H_INCLUDED
#define UTILS_INCLUDE_TRAJECTORY_H_INCLUDED

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
  SerializationMetadata()
    : protocolVersion(0),
      rank(0),
      trajectoryType("unspecified")
  {};
  int protocolVersion;
  int rank;
  std::string trajectoryType;

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & protocolVersion & rank & trajectoryType;}
#endif // DO_NOT_USE_BOOST_SERIALIZATION
};

namespace details
{

void writeNextArchive(std::ofstream*, const std::ostringstream&);
void readNextArchive(std::ifstream&, std::istringstream&);
SerializationMetadata readMeta(std::ifstream& ifs, bool reset);

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
  
  cpputils::oarchive& writeMeta(cpputils::oarchive& oar)  const {return writeMeta_v(oar);}
  
  virtual ~Trajectory() {}

private:
  virtual void            evolve_v(double)       = 0; // A step of exactly deltaT
  virtual double         getTime_v()       const = 0;
  virtual double        getDtDid_v()       const = 0;
  
  virtual std::ostream& displayParameters_v(std::ostream&         ) const = 0;
  virtual std::ostream& display_v          (std::ostream&, int    ) const = 0;
  virtual std::ostream& displayKey_v       (std::ostream&, size_t&) const = 0;

  virtual cpputils::iarchive&  readState_v(cpputils::iarchive&)       = 0;
  virtual cpputils::oarchive& writeState_v(cpputils::oarchive&) const = 0;
  virtual cpputils::oarchive& writeMeta_v (cpputils::oarchive& oar) const {SerializationMetadata meta; return oar & meta;}

  virtual std::ostream& logOnEnd_v(std::ostream& os) const {return os;}
  
};


///////////
//
// Adaptive
//
///////////


template<typename A>
class Adaptive : public virtual Trajectory 
{
public:
  // Some parameter-independent code could still be factored out, but probably very little
  
  typedef evolved::Evolved<A> Evolved;

  void step(double deltaT) {step_v(deltaT);}
  
  virtual ~Adaptive() {}

protected:
  Adaptive(A&, typename Evolved::Derivs, double, double, double    , const A&, const evolved::Maker<A>&);

  Adaptive(A&, typename Evolved::Derivs, double, const ParsEvolved&, const A&, const evolved::Maker<A>&);

  typedef typename Evolved::ConstPtr ConstPtr;
  typedef typename Evolved::Ptr Ptr;
  
  const ConstPtr getEvolved() const {return ConstPtr(evolved_);}
  const      Ptr getEvolved()       {return          evolved_ ;}

  double getDtTry() const {return evolved_->getDtTry();}

  std::ostream& displayParameters_v(std::ostream&) const;

  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       {return iar & *evolved_;}
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const {return oar & *evolved_;}

private:
  double getDtDid_v() const {return evolved_->getDtDid();}

  void evolve_v(double deltaT) {evolved::evolve<Adaptive>(*this,deltaT);}

  double getTime_v() const {return evolved_->getTime();}

  virtual void step_v(double deltaT) = 0;
  // Prefer purely virtual functions, so that there is no danger of forgetting to override them. Very few examples anyway for a trajectory wanting to perform only a step of Evolved.
  // (Only Simulated, but neither Master, nor MCWF_Trajectory)

  const Ptr evolved_;

};


} // trajectory


#endif // UTILS_INCLUDE_TRAJECTORY_H_INCLUDED
