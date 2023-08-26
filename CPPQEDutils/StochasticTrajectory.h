// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Random.h"
#include "Trajectory.h"

#include <vector>


namespace cppqedutils::trajectory {


// TODO: define concept of stochastic_trajectory
  

/// Aggregate of parameters pertaining to stochastic simulations
/** \copydetails ParsRun */
template <typename RandomEngine>
struct ParsStochastic : randomutils::Pars<RandomEngine,ode::Pars<>>
{
  size_t nTraj; ///< number of trajectories in case of ensemble averaging

  ParsStochastic(popl::OptionParser& op) : randomutils::Pars<RandomEngine,ode::Pars<>>{op} {
    add(op,"StochasticTrajectory",::parameters::_("nTraj","Number of trajectories",500,nTraj));
  }

};



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
    return std::ranges::fold_left_first (trajs | std::views::transform( [] (const SingleTrajectory& s) {return s.averaged();} ), std::plus{} ).value()/trajs.size();
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
  
  void advance(double deltaT, std::ostream& logStream) {for (auto& t : trajs_) cppqedutils::advance(t,deltaT,logStream);}

  /// An average of `dtDid`s from individual trajectories.
  double getDtDid() const
  {
    return std::ranges::fold_left_first( trajs_ | std::views::transform([] (const SingleTrajectory& t) {
      return t.getDtDid();
    }), std::plus{} ).value()/trajs_.size();
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
