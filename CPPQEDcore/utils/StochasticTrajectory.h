/// \briefFile{Defines classes related to stochastic evolution}
// -*- C++ -*-
#ifndef UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
#define UTILS_STOCHASTICTRAJECTORY_H_INCLUDED


#include "StochasticTrajectoryFwd.h"

#include "Trajectory.h"
#include "Randomized.h"


#include <boost/ptr_container/ptr_vector.hpp>


namespace trajectory {


/// The very general concept of an averageable trajectory
/**
 * Besides being a Trajectory, it can report a certain set of quantities, which are to be averaged either
 * - along a single Averageable trajectory (time average) or
 * - over several instances of Averageable all evolved to a certain time instant (ensemble average)
 * 
 * \tparam T the type condensing the quantities to be averaged. No implicit interface assumed @ this point.
 * Possible models: `double` or `complex` for a single c-number quantity; an \refStdCppConstruct{std::valarray,valarray/valarray}, or a quantumdata::DensityOperator
 * 
 * \todo implement general time averaging along the lines discussed in [this tracker](https://sourceforge.net/p/cppqed/feature-requests/1/#f3b3)
 * 
 */
template<typename T> 
class Averageable : public virtual Trajectory
{
public:
  virtual ~Averageable() {}

  const T toBeAveraged() const {return toBeAveraged_v();} ///< returns the set of quantities condensed in a variable of type `T` that are “to be averaged”

private:
  virtual const T toBeAveraged_v() const = 0;

};


/// Represents a trajectory that has both adaptive ODE evolution and noise
/** In the language of the framework, this means that the class simply connects Adaptive and Averageable while storing a randomized::Randomized instant for the convenience of derived classes. */
template<typename A, typename T>
class Stochastic : public Adaptive<A>, public Averageable<T>
{
public:
  virtual ~Stochastic() {}

private:
  typedef Adaptive<A> Base;

  typedef typename Base::Evolved Evolved;

  typedef randomized::Randomized::Ptr RandomizedPtr;

protected:
  /// \name Constructors
  //@{
  Stochastic(A&, typename Evolved::Derivs, double dtInit, 
             double epsRel, double epsAbs, const A& scaleAbs, 
             const evolved::Maker<A>&,
             unsigned long seed,
             bool noise,
             const randomized::Maker&); ///< Straightforward constructor combining the construction of Adaptive and randomized::Randomized

  Stochastic(A&, typename Evolved::Derivs, double dtInit,
             const A& scaleAbs, const ParsStochastic&,
             const evolved::Maker<A>&,
             const randomized::Maker&); ///< \overload
  //@}
  
  /// \name Getters
  //@{
  const RandomizedPtr getRandomized() const {return randomized_;}
  bool                isNoisy      () const {return isNoisy_   ;}
  //@}
  
  std::ostream& displayParameters_v(std::ostream&) const;
  
  /// \name Serialization
  //@{
  cpputils::iarchive&  readStateMore_v(cpputils::iarchive& iar)       {return iar & *randomized_;}
  cpputils::oarchive& writeStateMore_v(cpputils::oarchive& oar) const {return oar & *randomized_;}
  //@}
  
private:
  const unsigned long seed_ ;
  const bool          isNoisy_;

  const RandomizedPtr randomized_;

};


///////////
//
// Ensemble
//
///////////


namespace ensemble {

/// A base-class to Ensemble with customized behaviour according to the type of `T`
template<typename T>
class Base {};


/// Specialization for reference types
template<typename T>
class Base<T&>
{
public:
  typedef T& ToBeAveragedType;

  const ToBeAveragedType getInitializedToBeAveraged() const {return getInitializedToBeAveraged_v();}

  virtual ~Base() {}

private:
  virtual const ToBeAveragedType getInitializedToBeAveraged_v() const = 0;
  
};


} // ensemble

/*
  
  The implicit interface:
  
      reference case : T_ELEM must be addable to T via an addTo function.

  non-reference case : T must be constructible from a T_ELEM

*/


template<typename T, typename T_ELEM>
class Ensemble : public Averageable<T>, public ensemble::Base<T>
{
public:
  typedef Averageable<T     > Base;
  typedef Averageable<T_ELEM> Elem;

  typedef T ToBeAveragedType;

  const ToBeAveragedType averageInRange(size_t begin, size_t n) const;
  // Averages only in a range begin..begin+n-1. Earlier, this was called toBeAveraged, too, but it is not good to redefine an inherited non-virtual function.

  virtual ~Ensemble() {}

  typedef boost::ptr_vector<Elem> Impl;
  // We use a vector in order that individual trajectories are addressed more easily.

  const Impl& getTrajs() const {return trajs_;}

protected:
  typedef std::auto_ptr<Impl> Ptr;
  
  Ensemble(Ptr trajs, bool log) : trajs_(trajs), log_(log) {}

#define FOR_EACH_function(f) for_each(trajs_,bind(&Elem::f,_1,boost::ref(ios))); return ios;
  
  cpputils::iarchive&  readState_v(cpputils::iarchive& ios)       {FOR_EACH_function( readState)}
  cpputils::oarchive& writeState_v(cpputils::oarchive& ios) const {FOR_EACH_function(writeState)}

private:
  std::ostream& logOnEnd_v(std::ostream& ios) const {FOR_EACH_function(logOnEnd)}

#undef FOR_EACH_function
  
  void evolve_v(double deltaT);

  double getTime_v() const {return trajs_.front().getTime();}

  std::ostream& displayParameters_v(std::ostream&) const;

  double getDtDid_v() const;
  // An average of getDtDid()-s from individual trajectories.

  const ToBeAveragedType toBeAveraged_v() const {return averageInRange(0,trajs_.size());}

  Impl trajs_; // cannot be const because ptr_vector “propagates constness” (very correctly)

  const bool log_;

};



namespace ensemble {



// Implementation of the traits class for the most commonly used case:
template<typename T, typename T_ELEM>
class Traits
{
public:
  typedef Ensemble<T,T_ELEM> EnsembleType;

  typedef typename EnsembleType::Elem             Elem            ;
  typedef typename EnsembleType::Impl             Impl            ;
  typedef typename EnsembleType::ToBeAveragedType ToBeAveragedType;

  static const ToBeAveragedType averageInRange(typename Impl::const_iterator, typename Impl::const_iterator, const EnsembleType&);

};

} // ensemble


} // trajectory


#endif // UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
