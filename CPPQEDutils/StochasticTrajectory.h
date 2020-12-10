// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines classes related to stochastic evolution}
#ifndef CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED

#include "ParsStochasticTrajectory.h"
#include "Randomized.h"
#include "Trajectory.tcc"

#include "Conversions.h"

#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/numeric.hpp>

#include <boost/progress.hpp>

#include <boost/mpl/identity.hpp>



namespace mpl=boost::mpl;


namespace trajectory {


/// Contains important helpers for averaging (over time or ensemble)
namespace averaging {

/// Whereby a type T can safely and efficiently be passed around (for storage, function parameter/return, ect.)
/**
 * by default, the handle type is the same as the type itself (for types that can be passed around by value)
 */
template<typename T>
struct HandleType : mpl::identity<T> {};

} // averaging


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
template<typename SA, typename T> 
class Averageable : public Trajectory<SA>
{
public:
  virtual ~Averageable() {}

  typedef typename averaging::HandleType<T>::type AveragedHandle;
  
  const AveragedHandle averaged() const {return averaged_v();} ///< returns the set of quantities condensed in a variable of type `T` that are “to be averaged”

private:
  virtual const AveragedHandle averaged_v() const = 0;

};


/// Represents a trajectory that has both adaptive ODE evolution and noise
/** In the language of the framework, this means that the class simply connects Adaptive and Averageable while storing a randomized::Randomized instant for the convenience of derived classes. */
template<typename SA, typename A, typename T>
class Stochastic : public Adaptive<A,Averageable<SA,T>>
{
public:
  virtual ~Stochastic() {}

private:
  typedef Adaptive<A,Averageable<SA,T>> Base;

  typedef typename Base::Evolved Evolved;

  typedef randomized::Randomized::Ptr RandomizedPtr;

protected:
  /// \name Constructors
  //@{
  /// Straightforward constructor combining the construction of Adaptive and randomized::Randomized
  Stochastic(A& y, typename Evolved::Derivs derivs,
             double dtInit,
             int logLevel,
             double epsRel, double epsAbs, const A& scaleAbs,
             const evolved::Maker<A>& makerE,
             unsigned long seed,
             bool n,
             const randomized::Maker& makerR)
    : Base(y,derivs,dtInit,logLevel,epsRel,epsAbs,scaleAbs,makerE),
      seed_(seed), isNoisy_(n), randomized_(makerR(seed)) {}

  /// \overload
  Stochastic(A& y, typename Evolved::Derivs derivs,
             double dtInit,
             const A& scaleAbs,
             const ParsStochastic& p,
             const evolved::Maker<A>& makerE,
             const randomized::Maker& makerR)
    : Stochastic(y,derivs,dtInit,p.logLevel,p.epsRel,p.epsAbs,scaleAbs,makerE,p.seed,p.noise,makerR) {}
  //@}
  
  /// \name Getters
  //@{
  const RandomizedPtr getRandomized() const {return randomized_;}
  bool                isNoisy      () const {return isNoisy_   ;}
  //@}
  
  std::ostream& streamParameters_v(std::ostream& os) const override
  {
    return Base::streamParameters_v(os)<<"Stochastic Trajectory Parameters: seed="<<seed_<<std::endl<<(isNoisy_ ? "" : "No noise.\n");
  }
  
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
template<typename SA, typename T, typename T_ELEM=T>
class Ensemble : public Averageable<SA,T>
{
private:
  typedef Averageable<SA,T     > Base;
  
public:
  typedef Averageable<SA,T_ELEM> Elem;

  /// The storage of the element trajectories is through a \refBoost{pointer-vector,ptr_container/doc/ptr_vector.html}
  /** This correctly handles the elements’s eventual polymorphy and allows for random access to individual trajectories (cf. averageInRange()) */
  typedef boost::ptr_vector<Elem> Trajectories;
  
  typedef T Averaged;
  
  typedef typename averaging::HandleType<T>::type AveragedHandle;

  /// Averages only in a range `begin..begin+n-1`
  /** It could be called `averaged` as well, but it is not good to redefine an \link Averageable::averaged inherited non-virtual function\endlink. */
  const AveragedHandle averageInRange(size_t begin, size_t n) const;

  virtual ~Ensemble() {}

protected:
  // static helpers to constructor
  typedef std::unique_ptr<Trajectories> Ptr;
  
  /// Generic constructor
  Ensemble(Ptr&& trajs, ///< the sequence of elements owned by the `ptr_vector`
           bool displayProgress ///< when true, a \refBoost{display of progress,timer/doc/index.html} through the element trajectories will show up on `cerr` at each step of time evolution
          ) : trajs_(std::move(trajs)), displayProgress_(displayProgress) {}

  /// Getter
  const Trajectories& getTrajectories() const {return trajs_;}

#define FOR_EACH_function(f) for (auto& t : trajs_) t.f(ios); return ios;
  
  /// \name Serialization
  //@{
  cpputils::iarchive&  readState_v(cpputils::iarchive& ios)       final {FOR_EACH_function( readState)}
  cpputils::oarchive& writeState_v(cpputils::oarchive& ios) const final {FOR_EACH_function(writeState)}
  //@}
  
private:
  std::ostream& logOnEnd_v(std::ostream& ios) const override {FOR_EACH_function(logOnEnd)}

#undef FOR_EACH_function
  
  void evolve_v(double deltaT) final
  {
    using namespace boost;
    
    if (displayProgress_) {
      progress_display pd(trajs_.size(),std::cerr);
      for (auto i=trajs_.begin(); i!=trajs_.end(); (++i, ++pd)) i->evolve(deltaT);
    }
    else
      for_each(trajs_,bind(&Trajectory<SA>::evolve,_1,deltaT));
  }

  double getTime_v() const final {return trajs_.front().getTime();}

  std::ostream& streamParameters_v(std::ostream& os) const final {return trajs_.front().streamParameters( os<<"Ensemble of "<<trajs_.size()<<" trajectories."<<std::endl );}

  /// An average of getDtDid()-s from individual trajectories.
  double getDtDid_v() const final {return accumulate(trajs_,0.,[] (double init, const Trajectory<SA>& t) {return init+t.getDtDid();})/size2Double(trajs_.size());}

  const AveragedHandle averaged_v() const final {return averageInRange(0,trajs_.size());}

  Trajectories trajs_; // cannot be const because ptr_vector “propagates constness” (very correctly)

  const bool displayProgress_;

};


namespace averaging {


/// Governs how to average up several `T_ELEM` types into a `T` type in the most efficient way (which is usually not with the naive addition operator)
/**
 * \tparam T the averaged type of the ensemble
 * \tparam T_ELEM the averaged type of the underlying Averageable instances
 * 
 * A generic (naive) implementation is provided for the traits class right away, assuming that `T_ELEM` is additive and dividable by a double, and that it can be converted into a `T`.
 * 
 * \note The wrapper-class solution is necessary here as the function parameter types cannot be inferred due to heavy type-dependence
 * 
 */
template<typename SA, typename T, typename T_ELEM>
struct AverageTrajectoriesInRange
{
  /// Naive generic implementation
  typedef typename Ensemble<SA,T,T_ELEM>::Trajectories::const_iterator CI;
  static const typename HandleType<T>::type _(CI begin, CI end)
  {
    using namespace boost;
    return accumulate(++begin,end,begin->averaged(),[] (const auto& init, const typename Ensemble<T,T_ELEM>::Elem & e) {return init + e.averaged();})/size2Double(end-begin);
  }
  
};


} // averaging


} // trajectory


template<typename SA, typename T, typename T_ELEM>
auto
trajectory::Ensemble<SA,T,T_ELEM>::averageInRange(size_t begin, size_t n) const -> const AveragedHandle
{
  return averaging::AverageTrajectoriesInRange<SA,T,T_ELEM>::_(trajs_.begin()+begin,trajs_.begin()+(begin+n));
}



#endif // CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
