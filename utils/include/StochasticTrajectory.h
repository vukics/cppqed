// -*- C++ -*-
/*
  Besides being a Trajectory, StochasticTrajectory takes into account
  the possibility of some noise during Evolution, generated by a
  Randomized.

  In this case the need for taking EnsembleAverage of a certain set of
  quantities over several such Trajectories arises. This set should be
  condensed into the template parameter T of StochasticTrajectory,
  on which the requirement is that it support vector-space operations.

  Eg a double or complex for a single c-number quantity, a valarray, or
  a DensityOperator as in C++QED

  The design is such that EnsembleStochastic is recursive: an
  EnsembleStochastic can act as an element in a larger
  EnsembleStochastic.

*/
#ifndef _STOCHASTIC_TRAJECTORY_H
#define _STOCHASTIC_TRAJECTORY_H


#include "StochasticTrajectoryFwd.h"

#include "Trajectory.h"
#include "Randomized.h"


#include<boost/ptr_container/ptr_list.hpp>


namespace trajectory {


///////////////////////////
//
// StochasticTrajectoryBase
//
///////////////////////////


template<typename T> 
class StochasticTrajectoryBase : public virtual TrajectoryBase
{
public:
  virtual ~StochasticTrajectoryBase() {}

  virtual const T  toBeAveraged() const = 0;

};


///////////////////////
//
// StochasticTrajectory
//
///////////////////////

template<typename A, typename T>
class StochasticTrajectory : public Trajectory<A>, public StochasticTrajectoryBase<T>
{
public:
  typedef Trajectory<A> Base;

  typedef typename Base::Evolved Evolved;

  virtual void displayParameters() const;

  virtual ~StochasticTrajectory() {}

protected:
  typedef randomized::Randomized::SmartPtr RandomizedSmartPtr;

  StochasticTrajectory(A&, typename Evolved::Derivs, double dtInit, 
		       double epsRel, double epsAbs, const A& scaleAbs, 
		       const evolved::Maker<A>&,
		       unsigned long seed,
		       bool n,
		       const randomized::Maker&);

  StochasticTrajectory(A&, typename Evolved::Derivs, double dtInit,
		       const A& scaleAbs, const ParsStochasticTrajectory&,
		       const evolved::Maker<A>&,
		       const randomized::Maker&);

  const RandomizedSmartPtr getRandomized() const {return randomized_;}
  bool                     noise        () const {return isNoisy_   ;}

private:
  const unsigned long seed_ ;
  const bool          isNoisy_;

  RandomizedSmartPtr randomized_;

};


///////////////////////
//
// EnsembleTrajectories
//
///////////////////////


namespace details {

template<typename T>
class EnsembleBase {};


template<typename T>
class EnsembleBase<T&>
{
public:
  typedef T& TBA_Type;

  virtual TBA_Type getInitializedTBA() const = 0;

  virtual ~EnsembleBase() {}
  
};


}

/*
  
  The implicit interface:
  
      reference case : T_ELEM must be addable to T via an addTo function.

  non-reference case : T must be constructible from a T_ELEM

*/


template<typename T, typename T_ELEM>
class EnsembleTrajectories : public StochasticTrajectoryBase<T>, public details::EnsembleBase<T>
{
public:
  typedef StochasticTrajectoryBase<T     > Base;
  typedef StochasticTrajectoryBase<T_ELEM> Elem;

  typedef T TBA_Type;

  void evolve(double deltaT) const;

  double getTime() const {return trajs_.begin()->getTime();}

  const TBA_Type toBeAveraged() const;

  void displayParameters() const;

  virtual ~EnsembleTrajectories() {}

  typedef boost::ptr_list<Elem> Impl;

  const Impl& getTrajs() const {return trajs_;}

protected:
  typedef std::auto_ptr<Impl> SmartPtr;
  
  EnsembleTrajectories(SmartPtr trajs, bool log) : trajs_(trajs), log_(log) {}

private:
  double getDtDid() const;
  // An average of getDtDid()-s from individual trajectories.

  const Impl trajs_;

  const bool log_;

};



// Implementation of the traits class for the most commonly used case:

template<typename T, typename T_ELEM>
class EnsembleTrajectoriesTraits
{
public:
  typedef EnsembleTrajectories<T,T_ELEM> ET;

  typedef typename ET::Elem     Elem    ;
  typedef typename ET::Impl     Impl    ;
  typedef typename ET::TBA_Type TBA_Type;

  static const TBA_Type toBeAveraged(const ET& et);

};


} // trajectory

#include "impl/StochasticTrajectory.tcc"

#endif // _STOCHASTIC_TRAJECTORY_H
