// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Random.h"
#include "Trajectory.h"

#include <vector>


namespace cppqedutils::trajectory {


/// Aggregate of parameters pertaining to stochastic simulations
/** \copydetails ParsRun */
template <std::uniform_random_bit_generator RandomEngine>
struct ParsStochastic : randomutils::Pars<RandomEngine,ode::Pars<>>
{
  size_t nTraj; ///< number of trajectories in case of ensemble averaging

  ParsStochastic(popl::OptionParser& op) : randomutils::Pars<RandomEngine,ode::Pars<>>{op} {
    add(op,"StochasticTrajectory",::parameters::_("nTraj","Number of trajectories",500,nTraj));
  }

};


template <typename T>
concept stochastic = uniform_step<T> && requires (T&& t)
  {
    typename T::EnsembleAverageResult;
    typename T::EnsembleAverageElement;
    { averaged(t) } -> std::same_as<typename T::EnsembleAverageElement>;
  };



/// Governs how to average up several `SingleTrajectory` types in the most efficient way (which is sometimes not with the naive addition operator)
/**
 * \tparam SingleTrajectory type of a single trajectory in the ensemble
 * 
 * A generic naive implementation is provided for the traits class right away, assuming that `SingleTrajectory::EnsembleAverageElement` is additive and dividable by a double,
 * and that it can be converted into `SingleTrajectory::EnsembleAverageResult` via assignment.
 * 
 */
template<stochastic SingleTrajectory>
struct AverageTrajectories
{
  /// Naive generic implementation
  static const auto& _(typename SingleTrajectory::EnsembleAverageResult& , const std::vector<SingleTrajectory>& trajs)
  {
    return std::ranges::fold_left_first (trajs | std::views::transform( [] (const SingleTrajectory& s) {return averaged(s);} ), std::plus{} ).value()/trajs.size();
  }
  
};


template<stochastic SingleTrajectory>
struct InitializeEnsembleFromArrayOnlyArchive;


/// An ensemble of Averageable trajectories providing services for ensemble averaging and evolving the element trajectories serially
/**
 * \note Time averaging does not use stepsize-weighting, as experience has shown that this leads to worse convergence (similarly to quantumtrajectory::TimeAveragingMCWF_Trajectory).
 * \todo Stepsize-weighting could eventually be enabled as an option by a switch
 * 
 * The design is recursive: since Ensemble itself inherits from Averageable, it can act as an element in a larger Ensemble.
 * 
 * TODO: resurrect ensemble logging
 */
template<stochastic SingleTrajectory, typename TDP_Calculator /*, typename Logger*/>
class Ensemble
{
public:
  using TrajectoriesContainer = std::vector<SingleTrajectory>;
  
  using EnsembleAverageResult = typename SingleTrajectory::EnsembleAverageResult;
  using EnsembleAverageElement = EnsembleAverageResult; // for the recursive usage

  TrajectoriesContainer trajs;
  TDP_Calculator tdpCalculator;
//  Logger logger;
  mutable EnsembleAverageResult ensembleAverageResult;

  /// An average of `dtDid`s from individual trajectories.
  friend double getDtDid(const Ensemble& e)
  {
    return std::ranges::fold_left_first( e.trajs | std::views::transform([] (const SingleTrajectory& t) {
      return getDtDid(t);
    }), std::plus{} ).value()/e.trajs.size();
  }

  friend double getTime(const Ensemble& e) {return getTime(e.trajs.front());}
  
  // TODO: what to return here as log?
  friend auto advance(const Ensemble& e, double deltaT)
  {
    ::cppqedutils::LogTree res;
    for (auto& t : e.trajs) cppqedutils::advance(t,deltaT);
    return res;
  }

  friend auto logIntro(const Ensemble& e) { return ::cppqedutils::LogTree{{"Ensemble #traj",e.trajs.size()},{"Single-trajectory parameters",logIntro(e.trajs.front())}}; }

  friend auto logOutro(const Ensemble& e) { return logOutro(e.trajs.front()); }

  friend ::cppqedutils::LogTree dataStreamKey(const Ensemble& e) {return e.tdpCalculator.dataStreamKey();}

  friend auto temporalDataPoint(const Ensemble& e)
  {
    return e.tdpCalculator(getTime(e),averaged(e));
  }
  
  friend ::cppqedutils::iarchive& readFromArrayOnlyArchive(Ensemble& e, cppqedutils::iarchive& iar) {return InitializeEnsembleFromArrayOnlyArchive<SingleTrajectory>::_(e.trajs,iar);}
  
  template <typename Archive>
  friend auto& stateIO(Ensemble& e, Archive& ar)
  {
    for (auto& t : e.trajs) stateIO(t,ar); return ar;
  }
  
  /// This is what gets averaged when we have an Ensemble of Ensembles
  const EnsembleAverageResult& averaged(const Ensemble& e) const {return AverageTrajectories<SingleTrajectory>::_(e.ensembleAverageResult,e.trajs);}
  
};


} // cppqedutils::trajectory



template <typename SingleTrajectory, typename TDP_Calculator/*, typename Logger*/>
struct cppqedutils::trajectory::MakeSerializationMetadata<cppqedutils::trajectory::Ensemble<SingleTrajectory,TDP_Calculator/*,Logger*/>>
{
  static auto _()
  {
    auto metadataForSingle=MakeSerializationMetadata<SingleTrajectory>::_();
    return SerializationMetadata{metadataForSingle.typeID,"Ensemble of "+metadataForSingle.trajectoryID,metadataForSingle.rank};
  }
};
