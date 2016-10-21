// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED

#include "EnsembleMCWFFwd.h"

#include "DO_Display.h"
#include "MCWF_Trajectory.h"
#include "DensityOperator.h"

#include "Conversions.h"
#include "SmartPtr.h"


namespace quantumtrajectory {


namespace ensemble {

using namespace mcwf;

#define STATE_VECTORS(r) boost::ptr_vector<quantumdata::StateVector<r> >

#define BASE_class trajectory::Ensemble< quantumdata::DensityOperator<RANK>&, const quantumdata::StateVector<RANK>& >

/// Less templatized base for EnsembleMCWF \tparamRANK
template<int RANK>
class Base
  : private boost::base_from_member<STATE_VECTORS(RANK) >,
    public BASE_class
{
private:
  typedef STATE_VECTORS(RANK) StateVectors;

#undef  STATE_VECTORS

protected:
  typedef BASE_class Ensemble;

#undef  BASE_class

private:
  typedef typename Ensemble::Impl Trajectories;

  typedef boost::base_from_member<StateVectors> StateVectorsBase;

  typedef MCWF_Trajectory<RANK> Single;

public:
  typedef typename Single::StateVector    StateVector   ;
  typedef typename Single::StateVectorLow StateVectorLow; 

  typedef typename structure::QuantumSystem<RANK>::Ptr QuantumSystemPtr;

protected:
  /// Straightforward constructor
  Base(
       const StateVector& psi, ///< the (pure-state) initial condition
       QuantumSystemPtr sys, ///< the structure::QuantumSystem to be simulated
       const Pars& p, ///< parameters of the simulation (contains \link Pars::nTraj the number of trajectories\endlink)
       const StateVectorLow& =StateVectorLow()
       );

  const QuantumSystemPtr getQS() const {return qs_;}

private:
  virtual std::ostream& logOnEnd_v(std::ostream& os) const final;
  
  // static helpers to constructor
  // boost ptr_vector expects an auto_ptr in its interface, so we suppress the warning about auto_ptr being deprecated
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  static std::auto_ptr<StateVectors> stateVectors(const StateVector& psi, size_t nTraj);
  static std::auto_ptr<Trajectories> trajectories(StateVectors& psis, QuantumSystemPtr qs, const Pars& p, const StateVectorLow& scaleAbs);
#pragma GCC diagnostic pop

  
  virtual quantumdata::DensityOperator<RANK>& getInitializedDensityOperator_v() const final {rho_=0; return rho_;}

  mutable quantumdata::DensityOperator<RANK> rho_;

  const QuantumSystemPtr qs_;

  const size_t nBins_, nJumpsPerBin_;

};

} // ensemble


/// Derived from trajectory::Ensemble `<` quantumdata::DensityOperator `<RANK>& , const` quantumdata::StateVector `<RANK>& >`, it implements an ensemble of \link MCWF_Trajectory MCWF trajectories\endlink started from a pure-state initial condition
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
template<int RANK, typename V>
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
  /** \tparam SYS the physical system – can be any type convertible to structure::QuantumSystem::Ptr via cpputils::sharedPointerize */
  template<typename SYS>
  EnsembleMCWF(
               const StateVector& psi, ///< the (pure-state) initial condition used to initialize all the element \link MCWF_Trajectory MCWF trajectories\endlink
               const SYS& sys, ///< object representing the quantum system to be simulated
               const mcwf::Pars& p, ///< parameters of the simulation (contains \link mcwf::Pars::nTraj the number of trajectories\endlink)
               bool negativity, ///< governs whether entanglement should be calculated, cf. display_densityoperator::_, quantumdata::negPT
               const StateVectorLow& scaleAbs=StateVectorLow() ///< has the same role as `scaleAbs` in evolved::Maker::operator()
               )
    : Base(psi,cpputils::sharedPointerize(sys),p,scaleAbs), doDisplay_(structure::qsa<RANK>(this->getQS()),negativity) {}

private:
  virtual std::ostream& display_v   (std::ostream& os, int precision) const final {return doDisplay_.display   (this->getTime(),this->toBeAveraged(),os,precision);}
  virtual std::ostream& displayKey_v(std::ostream& os, size_t& i    ) const final {return doDisplay_.displayKey(os,i);}

  const DO_Display doDisplay_;

};


} // quantumtrajectory

/** \cond SPECIALIZATION */

namespace trajectory { namespace ensemble {


template<int RANK>
class Base<quantumdata::DensityOperator<RANK>&>
{
public:
  quantumdata::DensityOperator<RANK>& getInitializedDensityOperator() const {return getInitializedDensityOperator_v();}

private:
  virtual quantumdata::DensityOperator<RANK>& getInitializedDensityOperator_v() const = 0;
  
};


template<int RANK>
class Traits<quantumdata::DensityOperator<RANK>&, const quantumdata::StateVector<RANK>&>
{
public:
  typedef Ensemble<quantumdata::DensityOperator<RANK>&, const quantumdata::StateVector<RANK>&> EnsembleType;

  typedef typename EnsembleType::Elem             Elem            ;
  typedef typename EnsembleType::Impl             Impl            ;
  typedef typename EnsembleType::ToBeAveragedType ToBeAveragedType;

  /// assumes that quantumdata::StateVector features a member quantumdata::StateVector::addTo()
  static const ToBeAveragedType averageInRange(typename Impl::const_iterator begin, typename Impl::const_iterator end, const EnsembleType& et)
  {
    ToBeAveragedType res(et.getInitializedDensityOperator());
    
    for (auto i=begin; i!=end; i++) i->toBeAveraged().addTo(res);

    return res/=size2Double(end-begin);

  }


};

 
} } // trajectory::ensemble

/** \endcond */

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
