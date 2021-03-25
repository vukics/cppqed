// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines classes related to stochastic evolution}
#ifndef CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED

#include "Random.h"
#include "Trajectory.h"

#include "Conversions.h"

#include <boost/range/numeric.hpp>

#include <vector>


namespace cppqedutils::trajectory {


/// Aggregate of parameters pertaining to stochastic simulations
/** \copydetails ParsRun */
template <typename RandomEngine>
struct ParsStochastic : randomutils::Pars<RandomEngine,ode_engine::Pars<>>
{
  /// whether the noise should be on or off
  /**
   * (if it makes sense to turn it off at all for a concrete Stochastic
   * – e.g. for a \link quantumtrajectory::MCWF_Trajectory Monte Carlo wave-function trajectory\endlink, turning off the noise means simply to disable quantum jumps)
   */
  bool &noise;
  
  size_t &nTraj; ///< number of trajectories in case of ensemble averaging

  ParsStochastic(parameters::Table& p, const std::string& mod="")
    : randomutils::Pars<RandomEngine,ode_engine::Pars<>>{p,mod},
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
/*
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
*/

/// Represents a trajectory that has both adaptive ODE evolution and noise
/** Simply connects Adaptive and Averageable while storing a RandomEngine instant for the convenience of derived classes. *//*
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
  Stochastic(A&& y, typename Evolved::Derivs derivs,
             double dtInit,
             int logLevel,
             double epsRel, double epsAbs, const A& scaleAbs,
             const evolved::Maker<A>& makerE,
             bool isNoisy,
             RandomEngine&& re)
    : Base{std::forward<A>(y),derivs,dtInit,logLevel,epsRel,epsAbs,scaleAbs,makerE},
      isNoisy_(isNoisy),
      reInitStateDescriptor_{ [re{std::move(re)}]() {std::ostringstream oss; oss<<re; return oss.str();}() },
      re_(std::forward<RandomEngine>(re)) {}

  /// \overload
  Stochastic(A&& y, typename Evolved::Derivs derivs,
             double dtInit,
             const A& scaleAbs,
             const ParsStochastic<RandomEngine>& p,
             const evolved::Maker<A>& makerE)
    : Stochastic{std::forward<A>(y),derivs,dtInit,p.logLevel,p.epsRel,p.epsAbs,scaleAbs,makerE,p.noise,
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
  cppqedutils::iarchive&  readStateMore_v(cppqedutils::iarchive& iar)       override {return iar & re_;}
  cppqedutils::oarchive& writeStateMore_v(cppqedutils::oarchive& oar) const override {return oar & re_;}
  //@}
  
private:
  const bool isNoisy_;
  
  const std::string reInitStateDescriptor_;

  RandomEngine re_;

};
*/


/// Governs how to average up several `SingleTrajectory` types in the most efficient way (which is sometimes not with the naive addition operator)
/**
 * \tparam SingleTrajectory type of a single trajectory in the ensemble
 * 
 * A generic naive implementation is provided for the traits class right away, assuming that `SingleTrajectory::EnsembleAverageElement` is additive and dividable by a double,
 * and that it can be converted into `SingleTrajectory::EnsembleAverageResult` via assignment.
 * 
 */
template<typename SingleTrajectory>
struct AverageTrajectories
{
  /// Naive generic implementation
  static const auto& _(typename SingleTrajectory::EnsembleAverageResult& res, const std::vector<SingleTrajectory>& trajs)
  {
    using namespace boost;
    return res = accumulate(++trajs.begin(),trajs.end(),trajs.begin()->averaged(),[] (const auto& init, const SingleTrajectory& s) {
      return init + s.averaged();
    }
    )/size2Double(trajs.size());
  }
  
};


template<typename SingleTrajectory>
struct InitializeEnsembleFromArrayOnlyArchive;


/// An ensemble of Averageable trajectories providing services for ensemble averaging and evolving the element trajectories serially
/**
 * \note Time averaging does not use stepsize-weighting, as experience has shown that this leads to worse convergence (similarly to quantumtrajectory::TimeAveragingMCWF_Trajectory).
 * 
 * \todo Stepsize-weighting could eventually be enabled as an option by a switch
 * 
 * The design is recursive: since Ensemble itself inherits from Averageable, it can act as an element in a larger Ensemble.
 * 
 * \tparam ST Type of the elementary trajectories that form the Ensemble 
 * 
 * At the level of Ensemble, no implicit interface is assumed for `T` and `T_ELEM` since Ensemble treats variables of these types only via ensemble::Traits.
 * It is important that the way the averaged `T` will be calculated from the sequence of `T_ELEM`%s can be tailored
 * because it might happen that the application cannot afford to store temporaries of `T` (for such an example, cf. quantumtrajectory::EnsembleMCWF)
 * 
 */
template<typename ST, typename Streamer, typename Logger>
class Ensemble
{
public:
  using SingleTrajectory = ST;
  
  using TrajectoriesContainer = std::vector<SingleTrajectory>;
  
  using EnsembleAverageResult = typename SingleTrajectory::EnsembleAverageResult;
  
  using StreamedArray = typename SingleTrajectory::StreamedArray;
  
  Ensemble(Ensemble&&) = default; Ensemble& operator=(Ensemble&&) = default;

  /// Generic constructor
  template <typename TC, typename STR, typename LOG, typename ... Args>
  Ensemble(TC&& trajs, ///< the sequence of Singles (in general, there is no way to create different Singles from a given set of ctor parameters)
           STR&& streamer,
           LOG&& logger,
           Args&&... args)
  : trajs_{std::forward<TrajectoriesContainer>(trajs)},
    streamer_{std::forward<Streamer>(streamer)},
    logger_{std::forward<Logger>(logger)},
    ensembleAverageResult_{std::forward<Args>(args)...} {}

  auto getTime() const {return trajs_.front().getTime();}
  
  void evolve(double deltaT, std::ostream& logStream) {for (auto& t : trajs_) cppqedutils::evolve(t,deltaT,logStream);}

  /// An average of `dtDid`s from individual trajectories.
  double getDtDid() const {return boost::accumulate(trajs_,0.,[] (double init, const SingleTrajectory& t) {
      return init+t.getDtDid();
    })/size2Double(trajs_.size());
  }
  
  std::ostream& streamParameters(std::ostream& os) const {return trajs_.front().streamParameters( os<<"Ensemble of "<<trajs_.size()<<" trajectories."<<std::endl );}

  auto stream(std::ostream& os, int precision) const
  {
    return streamer_(getTime(),averaged(),os,precision);
  }
  
  auto& readFromArrayOnlyArchive(cppqedutils::iarchive& iar) {InitializeEnsembleFromArrayOnlyArchive<SingleTrajectory>::_(trajs_,iar); return iar;}
  
  template <typename Archive>
  auto& stateIO(Archive& ar) {for (auto& t : trajs_) t.stateIO(ar); return ar;}
  
  std::ostream& streamKey(std::ostream& os) const {size_t i=3; return streamer_.streamKey(os,i);}

  std::ostream& logOnEnd(std::ostream& os) const {logger_(trajs_,os); return os;};

  /// This is what gets averaged when we have an ensemble of Ensembles
  const EnsembleAverageResult& averaged() const {return AverageTrajectories<SingleTrajectory>::_(ensembleAverageResult_,trajs_);}
  
private:
  TrajectoriesContainer trajs_;
  
  Streamer streamer_;
  
  Logger logger_;
  
  mutable EnsembleAverageResult ensembleAverageResult_;

};


/// Deduction guide:
template <typename SingleTrajectory, typename Streamer, typename Logger, typename ... Args>
Ensemble(std::vector<SingleTrajectory>, Streamer, Logger, Args...) -> Ensemble<SingleTrajectory,Streamer,Logger>;


} // cppqedutils::trajectory



template <typename SingleTrajectory, typename Streamer, typename Logger>
struct cppqedutils::trajectory::MakeSerializationMetadata<cppqedutils::trajectory::Ensemble<SingleTrajectory,Streamer,Logger>>
{
  static auto _()
  {
    auto metadataForSingle=MakeSerializationMetadata<SingleTrajectory>::_();
    return SerializationMetadata{metadataForSingle.typeID,"Ensemble of "+metadataForSingle.trajectoryID,metadataForSingle.rank};
  }
};



#endif // CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_H_INCLUDED
