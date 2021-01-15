// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines classes related to stochastic evolution}
#ifndef CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED

#include "ParsTrajectory.h"
#include "Random.h"
#include "Trajectory.tcc"

#include "Conversions.h"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/numeric.hpp>



namespace trajectory {


/// Aggregate of parameters pertaining to stochastic simulations
/** \copydetails ParsRun */
template <typename RandomEngine>
struct ParsStochastic : randomutils::Pars<RandomEngine,ParsEvolved>
{
  /// whether the noise should be on or off
  /**
   * (if it makes sense to turn it off at all for a concrete Stochastic
   * – e.g. for a \link quantumtrajectory::MCWF_Trajectory Monte Carlo wave-function trajectory\endlink, turning off the noise means simply to disable quantum jumps)
   */
  bool &noise;
  
  size_t &nTraj; ///< number of trajectories in case of ensemble averaging

  ParsStochastic(parameters::Table& p, const std::string& mod="")
    : randomutils::Pars<RandomEngine,ParsEvolved>{p,mod},
      noise(p.add("noise",mod,"Switching noise on/off",true)),
      nTraj(p.add("nTraj",mod,"Number of trajectories",size_t(500)))
  {}
      
};



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
protected:
  Averageable() = default;
  Averageable(Averageable&&) = default; Averageable& operator=(Averageable&&) = default;

public:
  using ToBeAveraged=T;
  
  virtual ~Averageable() {}

  T averaged() const {return averaged_v();} ///< returns the set of quantities condensed in a variable of type `T` that are “to be averaged”

private:
  virtual T averaged_v() const = 0;

};


/// Represents a trajectory that has both adaptive ODE evolution and noise
/** Simply connects Adaptive and Averageable while storing a RandomEngine instant for the convenience of derived classes. */
template<typename SA, typename A, typename T, typename RandomEngine>
class Stochastic : public Adaptive<A,Averageable<SA,T>>
{
protected:
  Stochastic(Stochastic&&) = default; Stochastic& operator=(Stochastic&&) = default;

public:
  virtual ~Stochastic() {}

private:
  typedef Adaptive<A,Averageable<SA,T>> Base;

  typedef typename Base::Evolved Evolved;

protected:
  /// \name Constructors
  //@{
  /// Straightforward constructor combining the construction of Adaptive and a random engine
  Stochastic(A& y, typename Evolved::Derivs derivs,
             double dtInit,
             int logLevel,
             double epsRel, double epsAbs, const A& scaleAbs,
             const evolved::Maker<A>& makerE,
             bool isNoisy,
             RandomEngine&& re)
    : Base{y,derivs,dtInit,logLevel,epsRel,epsAbs,scaleAbs,makerE},
      isNoisy_(isNoisy),
      reInitStateDescriptor_{ [re{std::move(re)}]() {std::ostringstream oss; oss<<re; return oss.str();}() },
      re_(std::forward<RandomEngine>(re)) {}

  /// \overload
  Stochastic(A& y, typename Evolved::Derivs derivs,
             double dtInit,
             const A& scaleAbs,
             const ParsStochastic<RandomEngine>& p,
             const evolved::Maker<A>& makerE)
    : Stochastic{y,derivs,dtInit,p.logLevel,p.epsRel,p.epsAbs,scaleAbs,makerE,p.noise,
                 randomutils::streamOfOrdo(p)} {}
  //@}
  
  /// \name Getters
  //@{
public:
  auto& getRandomEngine() {return re_;}
  
protected:
  auto isNoisy() const {return isNoisy_;}
  //@}
  
  std::ostream& streamParameters_v(std::ostream& os) const override
  {
    return Base::streamParameters_v(os)<<"Stochastic trajectory random engine, initial state: "<<randomutils::EngineID_v<RandomEngine>
                                       <<" "<<reInitStateDescriptor_<<std::endl<<(isNoisy_ ? "" : "No noise.\n");
  }
  
  /// \name Serialization
  //@{
  cpputils::iarchive&  readStateMore_v(cpputils::iarchive& iar)       override {return iar & re_;}
  cpputils::oarchive& writeStateMore_v(cpputils::oarchive& oar) const override {return oar & re_;}
  //@}
  
private:
  const bool isNoisy_;
  
  const std::string reInitStateDescriptor_;

  RandomEngine re_;

};



/// An ensemble of Averageable trajectories providing services for ensemble averaging and evolving the element trajectories serially
/**
 * \note Time averaging does not use stepsize-weighting, as experience has shown that this leads to worse convergence (similarly to quantumtrajectory::TimeAveragingMCWF_Trajectory).
 * 
 * \todo Stepsize-weighting could eventually be enabled as an option by a switch
 * 
 * The design is recursive: since Ensemble itself inherits from Averageable, it can act as an element in a larger Ensemble.
 * 
 * \tparam ST Type of the elementary trajectories that form the Ensemble 
 * \tparam T The type condensing the quantities to be averaged for Ensemble in its function as an Averageable
 * 
 * At the level of Ensemble, no implicit interface is assumed for `T` and `T_ELEM` since Ensemble treats variables of these types only via ensemble::Traits.
 * It is important that the way the averaged `T` will be calculated from the sequence of `T_ELEM`%s can be tailored
 * because it might happen that the application cannot afford to store temporaries of `T` (for such an example, cf. quantumtrajectory::EnsembleMCWF)
 * 
 */
template<typename ST, typename SA, typename T>
class Ensemble : public Averageable<SA,T>
{
protected:
  Ensemble(Ensemble&&) = default; Ensemble& operator=(Ensemble&&) = default;

private:
  typedef Averageable<SA,T> Base;
  
public:
  using Single=ST;
  
  using SingleToBeAveraged=typename Single::ToBeAveraged;

  typedef std::vector<Single> Trajectories;
  
  /// Averages only in a range `begin..begin+n-1`
  /** It could be called `averaged` as well, but it is not good to redefine an \link Averageable::averaged inherited non-virtual function\endlink. */
  T averageInRange(size_t begin, size_t n) const;

  virtual ~Ensemble() {}

protected:
  /// Generic constructor
  Ensemble(Trajectories&& trajs ///< the sequence of Singles (in general, there is no way to create different Singles from a given set of ctor parameters)
          ) : trajs_(std::move(trajs)) {}

  /// Getter
  Trajectories& getTrajectories() {return trajs_;}  
  const Trajectories& getTrajectories() const {return trajs_;}

  /// \name Serialization
  //@{
  cpputils::iarchive& readState_v(cpputils::iarchive& ia) final {for (auto& t : trajs_) t.readState(ia); return ia;}
  cpputils::oarchive& writeState_v(cpputils::oarchive& oa) const final {for (auto& t : trajs_) t.writeState(oa); return oa;}
  //@}
  
private:
  std::ostream& logOnEnd_v(std::ostream& os) const override {for (auto& t : trajs_) t.logOnEnd(os); return os;}

  void evolve_v(double deltaT, std::ostream& logStream) final { for (auto& t : trajs_) t.evolve(deltaT,logStream); }

  double getTime_v() const final {return trajs_.front().getTime();} ///< all Singles should have the same time

  std::ostream& streamParameters_v(std::ostream& os) const final {return trajs_.front().streamParameters( os<<"Ensemble of "<<trajs_.size()<<" trajectories."<<std::endl );}

  /// An average of getDtDid()-s from individual trajectories.
  double getDtDid_v() const final {return boost::accumulate(trajs_,0.,[] (double init, const Trajectory<SA>& t) {
    return init+t.getDtDid();
  }
  )/size2Double(trajs_.size());}

  T averaged_v() const final {return averageInRange(0,trajs_.size());}

  Trajectories trajs_;

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
template<typename ST, typename SA, typename T>
struct AverageTrajectoriesInRange
{
  /// Naive generic implementation
  typedef typename Ensemble<ST,SA,T>::Trajectories::const_iterator CI;
  static auto _(CI begin, CI end)
  {
    using namespace boost;
    return accumulate(++begin,end,begin->averaged(),[] (const auto& init, const ST& s) {
      return init + s.averaged();
    }
    )/size2Double(end-begin);
  }
  
};


} // averaging



} // trajectory


template<typename ST, typename SA, typename T>
T trajectory::Ensemble<ST,SA,T>::averageInRange(size_t begin, size_t n) const
{
  return averaging::AverageTrajectoriesInRange<ST,SA,T>::_(trajs_.begin()+begin,trajs_.begin()+(begin+n));
}



#endif // CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
