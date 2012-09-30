// -*- C++ -*-

// AdaptiveTrajectory is more than Evolved only in that it takes into account
// the need for communicating towards the user from time to time during
// Evolution.

// This is manifested in the abstract function Display.

// Consider copying of Trajectories

#ifndef UTILS_INCLUDE_TRAJECTORY_H_INCLUDED
#define UTILS_INCLUDE_TRAJECTORY_H_INCLUDED

#include "TrajectoryFwd.h"

#include "Exception.h"
#include "Evolved.h"

#include <iostream>

#include <boost/utility.hpp>


namespace trajectory {


inline void runDt(Trajectory&, double time, double deltaT, bool displayInfo);

template<typename A>
inline void run  (AdaptiveTrajectory<A> &, double time, int          , bool displayInfo);

template<typename A>
inline void evolve(AdaptiveTrajectory<A>&, const ParsTrajectory&);


class StoppingCriterionReachedException : public cpputils::Exception {};


class OutfileOpeningException : public cpputils::TaggedException
{
public:
  OutfileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

};


/////////////////
//
// Trajectory
//
/////////////////


class Trajectory : private boost::noncopyable
{
public:
  std::ostream& getOstream  () const {return   ostream_;}
  int           getPrecision() const {return precision_;}

  void display   () const;
  void displayKey() const;

  virtual void   evolve(double deltaT) const = 0; // A step of exactly deltaT

  virtual double getTime()             const = 0;

  virtual double getDtDid()            const = 0;

  virtual void   displayParameters()   const = 0;
  
  static  double factor() {return 10.;}

  virtual ~Trajectory();

protected:
  Trajectory(      std::ostream&, int);
  Trajectory(const std::string &, int);
  Trajectory(const ParsTrajectory&   );

  Trajectory() : ostream_(std::cerr), precision_(0) {assert(false);} // A dummy constructor, which should never be called

private:
  virtual void   displayMore   () const = 0; // LOGICALLY const
  virtual size_t displayMoreKey() const = 0;

  std::ostream& ostream_;

  const int precision_;

};

/////////////
//
// AdaptiveTrajectory
//
/////////////

template<typename A>
class AdaptiveTrajectory : public virtual Trajectory 
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

  void evolve(double deltaT) const {evolved::evolve<const AdaptiveTrajectory>(*this,deltaT);}

  virtual void displayParameters() const;

  virtual ~AdaptiveTrajectory() {}

protected:
  AdaptiveTrajectory(A&, typename Evolved::Derivs, double, double, double, const A&,
	     const evolved::Maker<A>&);

  AdaptiveTrajectory(A&, typename Evolved::Derivs, double, const A&, const ParsTrajectory&,
	     const evolved::Maker<A>&);

  typename Evolved::Ptr getEvolved() const {return evolved_;}

  double getDtDid() const {return evolved_->getDtDid();}
  double getDtTry() const {return evolved_->getDtTry();}

private:
  typename Evolved::Ptr evolved_;

};

namespace details {

template<typename T, typename L, typename D>
void run(T& traj, L l, D d, void (*doRun)(T&,L,D), bool timestep, bool displayInfo);

void doRun(Trajectory&, long   nDt , double deltaT);
// Evolves the system on an AdaptiveTrajectory for nDt display intervals deltaT

void doRun(Trajectory&, double time, double deltaT);
// Evolves the system on an AdaptiveTrajectory up to time T and Displays in every deltaT

template<typename A> 
void doRun(AdaptiveTrajectory<A>& traj, double time, int dc);

} // details


} // trajectory


#endif // UTILS_INCLUDE_TRAJECTORY_H_INCLUDED
