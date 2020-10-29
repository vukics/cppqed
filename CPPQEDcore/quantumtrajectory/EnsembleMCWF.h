// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED

#include "DO_Display.h"
#include "MCWF_Trajectory.tcc"
#include "DensityOperator.h"
#include "ParsMCWF_Trajectory.h"

#include "Conversions.h"


namespace quantumtrajectory {


namespace ensemble {

using namespace mcwf;

#define BASE_class trajectory::Ensemble< structure::AveragedCommon::Averages, quantumdata::DensityOperator<RANK> , quantumdata::StateVector<RANK> >

/// Less templatized base for EnsembleMCWF \tparamRANK
template<int RANK>
class Base
  : public BASE_class
{
protected:
  typedef BASE_class Ensemble;

#undef  BASE_class

private:
  typedef typename Ensemble::Trajectories Trajectories;

  typedef MCWF_Trajectory<RANK> Single;

public:
  typedef typename Single::StateVector    StateVector   ;
  typedef typename Single::StateVectorLow StateVectorLow; 

  typedef typename structure::QuantumSystem<RANK>::Ptr QuantumSystemPtr;

protected:
  /// Straightforward constructor
  Base(std::shared_ptr<const StateVector> psi, ///< the (pure-state) initial condition
       QuantumSystemPtr qs, ///< the structure::QuantumSystem to be simulated
       const Pars& p, ///< parameters of the simulation (contains \link Pars::nTraj the number of trajectories\endlink)
       const StateVectorLow& scaleAbs=StateVectorLow())
    : Ensemble([psi,qs,&p,&scaleAbs] {
        Trajectories res;
        p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

        for (size_t i=0; i<p.nTraj; (++i, ++p.seed) ) 
          res.push_back(new MCWF_Trajectory<RANK>(std::make_shared<StateVector>(*psi),qs,p,scaleAbs));

        return res.release();
      } (),p.logLevel<0),
      qs_(qs), nBins_(p.nBins), nJumpsPerBin_(p.nJumpsPerBin)
    {}

  const QuantumSystemPtr getQS() const {return qs_;}

private:
  std::ostream& logOnEnd_v(std::ostream& os) const final
  {
    LoggerList loggerList;
    for (auto& i : this->getTrajectories())
      if (const auto traj=dynamic_cast<const MCWF_Trajectory<RANK>*>(&i))
        loggerList.push_back(traj->getLogger());
  
    return displayLog(os,loggerList,nBins_,nJumpsPerBin_);
  }
  
  const QuantumSystemPtr qs_;

  const size_t nBins_, nJumpsPerBin_;

};

} // ensemble


/// Derived from trajectory::Ensemble `<` quantumdata::DensityOperator `<RANK> , ` quantumdata::StateVector `<RANK> >`, it implements an ensemble of \link MCWF_Trajectory MCWF trajectories\endlink started from a pure-state initial condition
/**
 * The class overrides trajectory::Trajectory::display in such a way that at the time instant of display,
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
 * \tparam V has the same function as the template parameter `V` in display_densityoperator::_, which class is used here for deriving quantum averages to display from the assembled density operator
 * 
 * \todo An additional constructor could be added to initialize the ensemble by a full density operator, which could be appropriately sampled.
 * 
 */
template<int RANK, typename V=tmptools::V_Empty>
class EnsembleMCWF : public ensemble::Base<RANK>
{
private:
  typedef ensemble::Base<RANK> Base;

  typedef display_densityoperator::_<RANK,V> DO_Display;

public:
  typedef typename Base::StateVectorLow StateVectorLow; 

  typedef typename Base::Ensemble Ensemble;

  typedef typename Base      ::    StateVector     StateVector;
  typedef typename DO_Display::DensityOperator DensityOperator;

  /// Templated constructor with the same idea as Master::Master
  EnsembleMCWF(
               std::shared_ptr<const StateVector> psi, ///< the (pure-state) initial condition used to initialize all the element \link MCWF_Trajectory MCWF trajectories\endlink
               typename structure::QuantumSystem<RANK>::Ptr sys, ///< represents the quantum system to be simulated
               const mcwf::Pars& p, ///< parameters of the simulation (contains \link mcwf::Pars::nTraj the number of trajectories\endlink)
               bool negativity, ///< governs whether entanglement should be calculated, cf. display_densityoperator::_, quantumdata::negPT
               const StateVectorLow& scaleAbs=StateVectorLow() ///< has the same role as `scaleAbs` in evolved::Maker::operator()
               )
    : Base(psi,sys,p,scaleAbs), doDisplay_(structure::qsa<RANK>(this->getQS()),negativity) {}

private:
  typename Base::DisplayReturnType display_v(std::ostream& os, int precision) const final
  {
    const auto averages{*this->averaged()};
    return doDisplay_.display(this->getTime(),averages,os,precision);
  }
  
  std::ostream& displayKey_v(std::ostream& os, size_t& i) const final {return doDisplay_.displayKey(os,i);}

  const DO_Display doDisplay_;

};


} // quantumtrajectory

/** \cond SPECIALIZATION */

namespace trajectory { namespace averaging {


template<int RANK>
struct HandleType<quantumdata::DensityOperator<RANK> > : mpl::identity<std::shared_ptr<quantumdata::DensityOperator<RANK> > > {};


template<typename DA, int RANK>
struct AverageTrajectoriesInRange<DA,quantumdata::DensityOperator<RANK>,quantumdata::StateVector<RANK> >
{
  typedef typename Ensemble<DA,quantumdata::DensityOperator<RANK>,quantumdata::StateVector<RANK> >::Trajectories::const_iterator CI;
  
  static const auto _(CI begin, CI end)
  {
    auto res(std::make_shared<quantumdata::DensityOperator<RANK> >(*begin->averaged()));
      
    for (auto i=begin+1; i<end; i++) i->averaged()->addTo(*res);

    *res/=size2Double(end-begin);
    return res;

  }

};

 
} } // trajectory::averaging

/** \endcond */

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
