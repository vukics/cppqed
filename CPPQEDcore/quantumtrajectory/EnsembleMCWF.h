// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED

#include "StreamDensityOperator.h"
#include "MCWF_Trajectory.h"
#include "DensityOperator.h"

#include "Conversions.h"

#include <boost/range/combine.hpp>// #include <ranges>

namespace quantumtrajectory {


namespace ensemble {

using namespace mcwf;

#define BASE_class trajectory::Ensemble<MCWF_Trajectory<RANK,RandomEngine>,structure::AveragedCommon::Averages, quantumdata::DensityOperator<RANK> >

/// Less templatized base for EnsembleMCWF \tparamRANK
template<int RANK, typename RandomEngine>
class Base
  : public BASE_class
{
protected:
  typedef BASE_class Ensemble;

#undef  BASE_class

private:
  typedef typename Ensemble::Trajectories Trajectories;

public:
  typedef typename Ensemble::Single Single;

  typedef typename Single::StateVector    StateVector   ;
  typedef typename Single::StateVectorLow StateVectorLow; 

  typedef typename structure::QuantumSystem<RANK>::Ptr QuantumSystemPtr;

protected:
  /// Straightforward constructor
  Base(const StateVector& psi, ///< the (pure-state) initial condition
       QuantumSystemPtr qs, ///< the structure::QuantumSystem to be simulated
       const Pars<RandomEngine>& p, ///< parameters of the simulation (contains \link Pars::nTraj the number of trajectories\endlink)
       const StateVectorLow& scaleAbs=StateVectorLow())
    : Ensemble([psi,qs,&p,&scaleAbs] {
        Trajectories res;
        
        p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

        for (size_t i=0; i<p.nTraj; ++i) {
          res.push_back(new Single(StateVector(psi),qs,p,scaleAbs));// emplace_back(StateVector(psi),qs,p,scaleAbs);
          randomutils::incrementForNextStream(p);
        }
        
        return res.release();
      } () ),
      qs_(qs), nBins_(p.nBins), nJumpsPerBin_(p.nJumpsPerBin)
    {}

  const QuantumSystemPtr getQS() const {return qs_;}

private:
  std::ostream& logOnEnd_v(std::ostream& os) const final
  {
    LoggerList loggerList;
    for (auto& i : this->getTrajectories())
      if (const auto traj=dynamic_cast<const Single*>(&i))
        loggerList.push_back(traj->getLogger());
  
    return streamLog(os,loggerList,nBins_,nJumpsPerBin_);
  }
  
  const QuantumSystemPtr qs_;

  const size_t nBins_, nJumpsPerBin_;

};

} // ensemble


/// Derived from trajectory::Ensemble `<` quantumdata::DensityOperator `<RANK> , ` quantumdata::StateVector `<RANK> >`, it implements an ensemble of \link MCWF_Trajectory MCWF trajectories\endlink started from a pure-state initial condition
/**
 * The class overrides trajectory::Trajectory::stream in such a way that at the time instant of streaming,
 * the ensemble-averaged density operator of the system gets assembled from the stochastic state vectors of the element \link MCWF_Trajectory MCWF trajectories\endlink as
 * \f[\rho_{\text{ensemble}}(t)=\frac1{\text{\scriptsize number of trajectories}}\sum_{i\in\{\text{set of trajectories}\}}\ket{\Psi_i(t)}\bra{\Psi_i(t)}.\f]
 * 
 * \note This is done via the quantumdata::StateVector::addTo function, so that always at most *one single* full density operator needs to be stored in memory.
 * This makes that this class can be used for systems of larger dimensionality than Master, whose underlying ODE driver needs to store several (typically 6–7) density-operator instants.
 * 
 * The set of state vectors and the element \link MCWF_Trajectory MCWF trajectories\endlink are *owned* by the class.
 * 
 * \note The class obviously does not inherit from trajectory::Adaptive (a single adaptive timestep would in general result in different stepsizes for the element trajectories),
 * so that it can be used only in \link trajectory::Trajectory::run deltaT-mode\endlink.
 * 
 * \tparam RANK arity of the Hilbert space
 * \tparam V has the same function as the template parameter `V` in stream_densityoperator::_, which class is used here for deriving quantum averages to stream from the assembled density operator
 * 
 * \todo An additional constructor could be added to initialize the ensemble by a full density operator, which could be appropriately sampled.
 * 
 */
template<int RANK, typename RandomEngine, typename V=tmptools::V_Empty>
class EnsembleMCWF : public ensemble::Base<RANK,RandomEngine>
{
private:
  typedef ensemble::Base<RANK,RandomEngine> Base;

  typedef stream_densityoperator::_<RANK,V> DO_Stream;

public:
  typedef typename Base::StateVectorLow StateVectorLow; 

  typedef typename Base::Ensemble Ensemble;

  typedef typename Base      ::    StateVector     StateVector;
  typedef typename DO_Stream::DensityOperator DensityOperator;

  /// Templated constructor with the same idea as Master::Master
  EnsembleMCWF(const StateVector& psi, ///< the (pure-state) initial condition used to initialize all the element \link MCWF_Trajectory MCWF trajectories\endlink
               typename structure::QuantumSystem<RANK>::Ptr sys, ///< represents the quantum system to be simulated
               const mcwf::Pars<RandomEngine>& p, ///< parameters of the simulation (contains \link mcwf::Pars::nTraj the number of trajectories\endlink)
               bool negativity, ///< governs whether entanglement should be calculated, cf. stream_densityoperator::_, quantumdata::negPT
               const StateVectorLow& scaleAbs=StateVectorLow() ///< has the same role as `scaleAbs` in evolved::Maker::operator()
               )
    : Base(psi,sys,p,scaleAbs), dos_(structure::qsa<RANK>(this->getQS()),negativity) {}

private:
  typename Base::StreamReturnType stream_v(std::ostream& os, int precision) const final
  {
    return dos_.stream(this->getTime(),this->averaged(),os,precision);
  }
  
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const final {return dos_.streamKey(os,i);}

  const DO_Stream dos_;

};


} // quantumtrajectory

/** \cond SPECIALIZATION */

namespace trajectory { namespace averaging {


template<int RANK, typename RandomEngine>
struct AverageTrajectoriesInRange<quantumtrajectory::MCWF_Trajectory<RANK,RandomEngine>,
                                  structure::AveragedCommon::Averages,
                                  quantumdata::DensityOperator<RANK>>
{
  typedef typename Ensemble<quantumtrajectory::MCWF_Trajectory<RANK,RandomEngine>,
                            structure::AveragedCommon::Averages,
                            quantumdata::DensityOperator<RANK>>::Trajectories::const_iterator CI;
  
  static const auto _(CI begin, CI end)
  {
    quantumdata::DensityOperator<RANK> res(begin->averaged());
      
    for (auto i=begin+1; i<end; i++) i->averaged().addTo(res);

    res/=size2Double(end-begin);
    return res;

  }

};

 
} } // trajectory::averaging

/** \endcond */

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
