// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED

#include "MasterFwd.h"

#include "Structure.h"

#include "DO_Display.h"
#include "QuantumTrajectory.h"
#include "Types.h"

#include "Exception.h"
#include "SmartPtr.h"
#include "VectorFromMatrixSliceIterator.h"

#include <boost/function.hpp>



namespace quantumtrajectory {


/// Auxiliary tools to Master
namespace master {


typedef trajectory::ParsEvolved Pars;


/// Thrown if the system is not applicable in Master-equation evolution
/**
 * \see structure::ExactCommon::applicableInMaster, \ref masterequationlimitations
 */
struct SystemNotApplicable : cpputils::Exception {};


/// The actual working base of Master in the case when blitzplusplus::basi::Iterator is used for implementing multi-matrix multiplications \tparamRANK
template<int RANK>
class Base : public QuantumTrajectory<RANK, trajectory::Adaptive<typename quantumdata::Types<RANK>::DensityOperatorLow,trajectory::Trajectory> >
{
public:
  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Exact        <RANK> Exact        ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian  ;
  typedef structure::Liouvillean  <RANK> Liouvillean  ;
  typedef structure::Averaged     <RANK> Averaged     ;

  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;
  typedef typename quantumdata::Types<RANK>::    StateVectorLow     StateVectorLow;

  typedef trajectory::Adaptive<DensityOperatorLow,trajectory::Trajectory> Adaptive;

  typedef quantumtrajectory::QuantumTrajectory<RANK, Adaptive> QuantumTrajectory;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  /// The actual function calculating the time derivative for \link evolved::Evolved ODE evolution\endlink
  void derivs(double, const DensityOperatorLow&, DensityOperatorLow&) const;

protected:
  using QuantumTrajectory::getQSW;

  Base(DensityOperator&, typename QuantumSystem::Ptr, const Pars&, const DensityOperatorLow& =DensityOperatorLow());

  typedef boost::function<void(                       StateVectorLow&)>  UnaryFunction;
  typedef boost::function<void(const StateVectorLow&, StateVectorLow&)> BinaryFunction;

  DensityOperator& rho_;

  const typename Averaged::Ptr getAv() const {return getQSW().getAv();}

private:
  void step_v(double) final;

  std::ostream& displayParameters_v(std::ostream&) const override;

  virtual void  unaryIter(                           DensityOperatorLow&,  UnaryFunction) const;
  virtual void binaryIter(const DensityOperatorLow&, DensityOperatorLow&, BinaryFunction) const;

  virtual const std::string addToParameterDisplay() const {return "";}

};



/// The actual working base of Master in the case when blitzplusplus::basi_fast::Iterator is used for implementing multi-matrix multiplications \tparamRANK
template<int RANK>
class BaseFast : public Base<RANK>
{
public:
  typedef typename Base<RANK>::QuantumSystem QuantumSystem;

  typedef typename Base<RANK>::DensityOperator DensityOperator;

  typedef typename Base<RANK>::DensityOperatorLow DensityOperatorLow;

  BaseFast(DensityOperator& rho, typename QuantumSystem::Ptr sys, const Pars& p, const DensityOperatorLow& scaleAbs=DensityOperatorLow())
    : Base<RANK>(rho,sys,p,scaleAbs), slicesData_(rho.getArray()) {}

private:
  typedef typename Base<RANK>:: UnaryFunction  UnaryFunction;
  typedef typename Base<RANK>::BinaryFunction BinaryFunction;

  void  unaryIter(                           DensityOperatorLow&,  UnaryFunction) const final;
  void binaryIter(const DensityOperatorLow&, DensityOperatorLow&, BinaryFunction) const final;

  const std::string addToParameterDisplay() const final {return " Fast Iteration.";}

  const blitzplusplus::SlicesData<2*RANK,blitzplusplus::vfmsi::LeftRight<RANK,blitzplusplus::vfmsi::Left> > slicesData_;

};



} // master



#define BASE_class boost::mpl::if_c<IS_FAST,master::BaseFast<RANK>,master::Base<RANK> >::type

/// An \link trajectory::Adaptive Adaptive\endlink trajectory class representing Master equation evolution from a \link quantumdata::DensityOperator density-operator\endlink initial condition
/**
 * \see \ref masterequation.

 * \note The ODE driver underlying this class needs to store several (typically 6–7, this depends on the chosen driver) density-operator instants.
 *
 * \tparam RANK arity of the Hilbert space
 * \tparam V has the same function as the template parameter `V` in display_densityoperator::_, which class is used here for deriving quantum averages to display from the evolved density operator
 * \tparam IS_FAST the class will use either blitzplusplus::basi::Iterator or blitzplusplus::basi_fast::Iterator to perform the multi-matrix multiplications, depending on this template argument
 * 
 */
template<int RANK, typename V, bool IS_FAST>
class Master : public BASE_class
{
public:
  typedef typename BASE_class Base;

#undef  BASE_class

  typedef typename Base::QuantumSystem QuantumSystem;

  typedef typename Base::DensityOperatorLow DensityOperatorLow; 

  typedef typename Base::DensityOperator DensityOperator;

  /// Templated constructor
  /** \tparam SYS the physical system – can be any type convertible to structure::QuantumSystem::Ptr via cpputils::sharedPointerize */
  template<typename SYS>
  Master(DensityOperator& rho, ///< the density operator to be evolved
         const SYS& sys, ///< object representing the quantum system
         const master::Pars& pt, ///< parameters of the evolution
         bool negativity, ///< governs whether entanglement should be calculated, cf. display_densityoperator::_, quantumdata::negPT
         const DensityOperatorLow& scaleAbs=DensityOperatorLow() ///< has the same role as `scaleAbs` in evolved::Maker::operator()
        )
    : Base(rho,cpputils::sharedPointerize(sys),pt,scaleAbs), doDisplay_(this->getAv(),negativity)
  {}
  
  const DensityOperator& getRho() const {return this->rho_;}

private:
  std::ostream& display_v   (std::ostream& os, int precision) const override {return doDisplay_.display   (this->getTime(),this->rho_,os,precision);}
  std::ostream& displayKey_v(std::ostream& os, size_t& i    ) const override {return doDisplay_.displayKey(os,i);}

  const std::string trajectoryID_v() const override {return "Master";}
  
  const display_densityoperator::_<RANK,V> doDisplay_;

};



} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
