// -*- C++ -*-

// Trajectory is more than Evolved only in that it takes into account
// the need for communicating towards the user from time to time during
// Evolution.

// This is manifested in the abstract function Display.

// Consider copying of Trajectories

#ifndef _TRAJECTORY_H
#define _TRAJECTORY_H

#include "TrajectoryFwd.h"

#include "Exception.h"
#include "Evolved.h"

#include<iostream>

#include<boost/utility.hpp>


namespace trajectory {


inline void runDt(TrajectoryBase&, double time, double deltaT, bool displayInfo);

template<typename A>
inline void run  (Trajectory<A> &, double time, int          , bool displayInfo);

template<typename A>
inline void evolve(Trajectory<A>&, const ParsTrajectory&);


class StoppingCriterionReachedException : public cpputils::Exception {};


class OutfileOpeningException : public cpputils::TaggedException
{
public:
  OutfileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

};


/////////////////
//
// TrajectoryBase
//
/////////////////

class TrajectoryBase : private boost::noncopyable
{
public:
  std::ostream& getOstream() const {return ostream_;}

  void display   () const;
  void displayKey() const;

  virtual void   evolve(double deltaT) const = 0; // A step of exactly deltaT

  virtual double getTime()             const = 0;

  virtual double getDtDid()            const = 0;

  virtual void   displayParameters()   const = 0;
  
  static  double factor() {return 10.;}

  virtual ~TrajectoryBase();

protected:
  TrajectoryBase(      std::ostream&, int);
  TrajectoryBase(const std::string &, int);
  TrajectoryBase(const ParsTrajectory&   );

  TrajectoryBase() : ostream_(std::cerr), precision_(0) {assert(false);} // A dummy constructor, which should never be called

private:
  virtual void   displayMore   (int)   const = 0; // LOGICALLY const
  virtual size_t displayMoreKey()      const = 0;

  std::ostream& ostream_;

  const int precision_;

};

/////////////
//
// Trajectory
//
/////////////

template<typename A>
class Trajectory : public virtual TrajectoryBase 
{
public:
  // Some parameter-independent code could still be factored out, but probably very little
  typedef evolved::Evolved<A> Evolved;

  double getTime() const {return evolved_->getTime();}

  // Preference for purely virtual functions, so that there is no
  // danger of forgetting to override them. Very few examples anyway
  // for a trajectory wanting to perform only a step of Evolved. (Only
  // Simulated, but neither Master, nor MCWF_Trajectory)
  virtual void step(double deltaT) const = 0;

  void evolve(double deltaT) const {evolved::evolve<const Trajectory>(*this,deltaT);}

  virtual void displayParameters() const;

  virtual ~Trajectory() {}

protected:
  Trajectory(A&, typename Evolved::Derivs, double, double, double, const A&,
	     const evolved::Maker<A>&);

  Trajectory(A&, typename Evolved::Derivs, double, const A&, const ParsTrajectory&,
	     const evolved::Maker<A>&);

  typename Evolved::SmartPtr getEvolved() const {return evolved_;}

  double getDtDid() const {return evolved_->getDtDid();}
  double getDtTry() const {return evolved_->getDtTry();}

private:
  typename Evolved::SmartPtr evolved_;

};



} // trajectory

#include "impl/Trajectory.tcc"

#endif // _TRAJECTORY_H
