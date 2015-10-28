// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines classes related to stochastic evolution}
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED


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
  
  std::ostream& displayParameters_v(std::ostream&) const override;
  
  /// \name Serialization
  //@{
  cpputils::iarchive&  readStateMore_v(cpputils::iarchive& iar)       override {return iar & *randomized_;}
  cpputils::oarchive& writeStateMore_v(cpputils::oarchive& oar) const override {return oar & *randomized_;}
  //@}
  
private:
  const unsigned long seed_ ;
  const bool          isNoisy_;

  const RandomizedPtr randomized_;

};



/// Contains important helpers to Ensemble
namespace ensemble {

/// A base-class to Ensemble with customized behaviour according to the type of `T`
template<typename T>
class Base {};

} // ensemble



/// An ensemble of Averageable trajectories providing services for ensemble averaging and evolving the element trajectories serially
/**
 * \note Time averaging does not use stepsize-weighting, as experience has shown that this leads to worse convergence (similarly to quantumtrajectory::TimeAveragingMCWF_Trajectory).
 * 
 * \todo Stepsize-weighting could eventually be enabled as an option by a switch
 * 
 * The design is recursive: since Ensemble itself inherits from Averageable, it can act as an element in a larger Ensemble.
 * 
 * The elements do not need to have the same type, they only need to have a common Averageable type as a base.
 * 
 * \tparam T The type condensing the quantities to be averaged for Ensemble in its function as an Averageable
 * \tparam T_ELEM Same for the Averageable%s 
 * 
 * At the level of Ensemble, no implicit interface is assumed for `T` and `T_ELEM` since Ensemble treats variables of these types only via ensemble::Traits.
 * It is important that the way the averaged `T` will be calculated from the sequence of `T_ELEM`%s can be tailored
 * because it might happen that the application cannot afford to store temporaries of `T` (for such an example, cf. quantumtrajectory::EnsembleMCWF)
 * 
 */
template<typename T, typename T_ELEM>
class Ensemble : public Averageable<T>, public ensemble::Base<T>
{
private:
  typedef Averageable<T     > Base;
  
public:
  typedef Averageable<T_ELEM> Elem;

  /// The storage of the element trajectories is through a \refBoost{pointer-vector,ptr_container/doc/ptr_vector.html}
  /** This correctly handles the elements’s eventual polymorphy and allows for random access to individual trajectories (cf. averageInRange()) */
  typedef boost::ptr_vector<Elem> Impl;
  
  typedef T ToBeAveragedType;

  /// Averages only in a range `begin..begin+n-1`
  /** It could be called `toBeAveraged` as well, but it is not good to redefine an \link Averageable::toBeAveraged inherited non-virtual function\endlink. */
  const ToBeAveragedType averageInRange(size_t begin, size_t n) const;

  virtual ~Ensemble() {}

protected:
  // static helpers to constructor
  // boost ptr_vector expects an auto_ptr in its interface, so we suppress the warning about auto_ptr being deprecated
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  typedef std::auto_ptr<Impl> Ptr;
#pragma GCC diagnostic pop
  
  /// Generic constructor
  Ensemble(Ptr trajs, ///< the sequence of elements owned by the `ptr_vector`
           bool displayProgress ///< when true, a \refBoost{display of progress,timer/doc/index.html} through the element trajectories will show up on `cerr` at each step of time evolution
          ) : trajs_(trajs), displayProgress_(displayProgress) {}

  /// Getter
  const Impl& getTrajectories() const {return trajs_;}

#define FOR_EACH_function(f) for_each(trajs_,bind(&Elem::f,_1,boost::ref(ios))); return ios;
  
  /// \name Serialization
  //@{
  virtual cpputils::iarchive&  readState_v(cpputils::iarchive& ios)       final {FOR_EACH_function( readState)}
  virtual cpputils::oarchive& writeState_v(cpputils::oarchive& ios) const final {FOR_EACH_function(writeState)}
  //@}
  
private:
  std::ostream& logOnEnd_v(std::ostream& ios) const override {FOR_EACH_function(logOnEnd)}

#undef FOR_EACH_function
  
  virtual void evolve_v(double deltaT) final;

  virtual double getTime_v() const final {return trajs_.front().getTime();}

  virtual std::ostream& displayParameters_v(std::ostream&) const final;

  virtual double getDtDid_v() const final;
  // An average of getDtDid()-s from individual trajectories.

  virtual const ToBeAveragedType toBeAveraged_v() const final {return averageInRange(0,trajs_.size());}

  Impl trajs_; // cannot be const because ptr_vector “propagates constness” (very correctly)

  const bool displayProgress_;

};


namespace ensemble {

/// %Traits class governing how to average up several `T_ELEM` types into a `T` type in the most efficient way (which is usually not with the naive addition operator)
/**
 * \tparam T the to-be-averaged type of the ensemble
 * \tparam T_ELEM the to-be-averaged type of the underlying Averageable instances
 * 
 * A generic (naive) implementation is provided for the traits class right away.
 * It assumes that `T_ELEM` is additive and dividable by a double, and that it can be converted into a `T`.
 * 
 */
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


#endif // CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
