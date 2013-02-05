// -*- C++ -*-
// Adaptive is more than Evolved only in that it takes into account the need for communicating towards the user from time to time during Evolution.
// This is manifested in the abstract function Display.
// Consider copying (cloning) of Trajectories

#ifndef UTILS_INCLUDE_TRAJECTORY_H_INCLUDED
#define UTILS_INCLUDE_TRAJECTORY_H_INCLUDED

#include "TrajectoryFwd.h"

#include "Exception.h"
#include "Evolved.h"

#include <boost/utility.hpp>

#ifndef   DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION


#include <iosfwd>



namespace trajectory {


void run(Trajectory &, double time, double deltaT, bool displayInfo);

void run(Trajectory &, long   nDt , double deltaT, bool displayInfo);

template<typename A>
void run(Adaptive<A>&, double time, int dc       , bool displayInfo);

void run(Trajectory &, const ParsRun&);

template<typename A>
void run(Adaptive<A>&, const ParsRun&);


class StoppingCriterionReachedException : public cpputils::Exception {};


class OutfileOpeningException : public cpputils::TaggedException
{
public:
  OutfileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

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

  void evolve(double deltaT) const {evolve_v(deltaT);} // A step of exactly deltaT

  double getTime() const {return getTime_v();}

  double getDtDid() const {return getDtDid_v();}

  std::ostream& displayParameters(std::ostream& os) const {return displayParameters_v(os);}
  
#ifndef   DO_NOT_USE_BOOST_SERIALIZATION
  typedef boost::archive::binary_iarchive iarchive;
  typedef boost::archive::binary_oarchive oarchive;
  
  iarchive&  readState(iarchive& iar)       {return  readState_v(iar);}
  oarchive& writeState(oarchive& oar) const {return writeState_v(oar);};
#endif // DO_NOT_USE_BOOST_SERIALIZATION
  
  virtual ~Trajectory() {}

private:
  virtual void            evolve_v(double) const = 0; // A step of exactly deltaT
  virtual double         getTime_v()       const = 0;
  virtual double        getDtDid_v()       const = 0;
  
  virtual std::ostream& displayParameters_v(std::ostream&) const = 0;

  virtual std::ostream& display_v          (std::ostream&, int    ) const = 0;
  virtual std::ostream& displayKey_v       (std::ostream&, size_t&) const = 0;

#ifndef   DO_NOT_USE_BOOST_SERIALIZATION
  virtual iarchive&  readState_v(iarchive&)       = 0;
  virtual oarchive& writeState_v(oarchive&) const = 0;
#endif // DO_NOT_USE_BOOST_SERIALIZATION

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

  void step(double deltaT) const {step_v(deltaT);}
  
  virtual ~Adaptive() {}

protected:
  Adaptive(A&, typename Evolved::Derivs, double, double, double    , const A&, const evolved::Maker<A>&);

  Adaptive(A&, typename Evolved::Derivs, double, const ParsEvolved&, const A&, const evolved::Maker<A>&);

  typename Evolved::Ptr getEvolved() const {return evolved_;}

  double getDtTry() const {return evolved_->getDtTry();}

  std::ostream& displayParameters_v(std::ostream&) const;

#ifndef   DO_NOT_USE_BOOST_SERIALIZATION
  iarchive&  readState_v(iarchive& iar)       {return iar & *evolved_;}
  oarchive& writeState_v(oarchive& oar) const {return oar & *evolved_;}
#endif // DO_NOT_USE_BOOST_SERIALIZATION

private:
  double getDtDid_v() const {return evolved_->getDtDid();}

  void evolve_v(double deltaT) const {evolved::evolve<const Adaptive>(*this,deltaT);}

  double getTime_v() const {return evolved_->getTime();}

  virtual void step_v(double deltaT) const = 0;
  // Prefer purely virtual functions, so that there is no danger of forgetting to override them. Very few examples anyway for a trajectory wanting to perform only a step of Evolved. (Only Simulated, but neither Master, nor MCWF_Trajectory)

  typename Evolved::Ptr evolved_;

};


} // trajectory


#endif // UTILS_INCLUDE_TRAJECTORY_H_INCLUDED
