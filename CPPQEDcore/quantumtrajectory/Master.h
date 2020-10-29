// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED

#include "Structure.h"

#include "DO_Display.h"
#include "QuantumTrajectory.h"
#include "Types.h"

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
struct SystemNotApplicable : std::runtime_error {SystemNotApplicable() : std::runtime_error("") {}};


/// The actual working base of Master in the case when blitzplusplus::basi::Iterator is used for implementing multi-matrix multiplications \tparamRANK
template<int RANK>
class Base : public QuantumTrajectory<RANK,
                                      trajectory::Adaptive<typename quantumdata::Types<RANK>::DensityOperatorLow,trajectory::Trajectory<typename structure::AveragedCommon::Averages>
                                                          > 
                                     >
{
public:
  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Exact        <RANK> Exact        ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian  ;
  typedef structure::Liouvillean  <RANK> Liouvillean  ;
  typedef structure::Averaged     <RANK> Averaged     ;

  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;
  typedef typename quantumdata::Types<RANK>::    StateVectorLow     StateVectorLow;

  typedef trajectory::Adaptive<DensityOperatorLow,trajectory::Trajectory<typename Averaged::Averages>> Adaptive;

  typedef quantumtrajectory::QuantumTrajectory<RANK, Adaptive> QuantumTrajectory;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;
  typedef std::shared_ptr<DensityOperator> DO_Ptr;

  /// The actual function calculating the time derivative for \link evolved::Evolved ODE evolution\endlink
  void derivs(double, const DensityOperatorLow&, DensityOperatorLow&) const;

protected:
  Base(DO_Ptr, typename QuantumSystem::Ptr, const Pars&, const DensityOperatorLow& =DensityOperatorLow());

  typedef boost::function<void(                       StateVectorLow&)>  UnaryFunction;
  typedef boost::function<void(const StateVectorLow&, StateVectorLow&)> BinaryFunction;

  const DO_Ptr rho_;

private:
  void step_v(double) final;

  std::ostream& displayParameters_v(std::ostream&) const override;

  virtual void  unaryIter(                           DensityOperatorLow&,  UnaryFunction) const;
  virtual void binaryIter(const DensityOperatorLow&, DensityOperatorLow&, BinaryFunction) const;

  virtual const std::string addToParameterDisplay() const {return "";}

};



} // master



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
template<int RANK, typename V=tmptools::V_Empty>
class Master : public master::Base<RANK>
{
public:
  typedef master::Base<RANK> Base;

  typedef typename Base::QuantumSystem QuantumSystem;

  typedef typename Base::DensityOperatorLow DensityOperatorLow; 

  typedef typename Base::DensityOperator DensityOperator;
  typedef std::shared_ptr<DensityOperator> DO_Ptr;

  /// Templated constructor
  Master(DO_Ptr rho, ///< the density operator to be evolved
         typename QuantumSystem::Ptr sys, ///< object representing the quantum system
         const master::Pars& pt, ///< parameters of the evolution
         bool negativity, ///< governs whether entanglement should be calculated, cf. display_densityoperator::_, quantumdata::negPT
         const DensityOperatorLow& scaleAbs=DensityOperatorLow() ///< has the same role as `scaleAbs` in evolved::Maker::operator()
        )
    : Base(rho,sys,pt,scaleAbs), doDisplay_(this->getAv(),negativity)
  {}
  
  const DO_Ptr getRho() const {return this->rho_;}

private:
  typename Base::DisplayReturnType display_v(std::ostream& os, int precision) const override {return doDisplay_.display(this->getTime(),*this->rho_,os,precision);}
  
  std::ostream& displayKey_v(std::ostream& os, size_t& i) const override {return doDisplay_.displayKey(os,i);}

  const std::string trajectoryID_v() const override {return "Master";}
  
  const display_densityoperator::_<RANK,V> doDisplay_;

};



} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
