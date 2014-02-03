/// \briefFileDefault
// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_MASTER_H_INCLUDED
#define QUANTUMTRAJECTORY_MASTER_H_INCLUDED

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


namespace master {


typedef trajectory::ParsEvolved Pars;


struct SystemNotApplicable  : cpputils::Exception {};
struct NoLiouvillean        : cpputils::Exception {};


template<int RANK>
class Base : public QuantumTrajectory<RANK, trajectory::Adaptive<typename quantumdata::Types<RANK>::DensityOperatorLow> >
{
public:
  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Exact        <RANK> Exact        ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian  ;
  typedef structure::Liouvillean  <RANK> Liouvillean  ;
  typedef structure::Averaged     <RANK> Averaged   ;

  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;
  typedef typename quantumdata::Types<RANK>::    StateVectorLow     StateVectorLow;

  typedef trajectory::Adaptive<DensityOperatorLow> Adaptive;

  typedef quantumtrajectory::QuantumTrajectory<RANK, Adaptive> QuantumTrajectory;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;
  
  using Adaptive::getEvolved; using Adaptive::getDtDid; using Adaptive::getTime;

  Base(DensityOperator&, typename QuantumSystem::Ptr, const Pars&, const DensityOperatorLow& =DensityOperatorLow());

  void derivs(double, const DensityOperatorLow&, DensityOperatorLow&) const;

protected:
  typedef boost::function<void(                       StateVectorLow&)>  UnaryFunction;
  typedef boost::function<void(const StateVectorLow&, StateVectorLow&)> BinaryFunction;

  DensityOperator& rho_;

  using QuantumTrajectory::getQSW;

  const typename Averaged::Ptr getAv() const {return getQSW().getAv();}

private:
  void              step_v(double);

  std::ostream& displayParameters_v(std::ostream&) const;

  virtual void  unaryIter(                           DensityOperatorLow&,  UnaryFunction) const;
  virtual void binaryIter(const DensityOperatorLow&, DensityOperatorLow&, BinaryFunction) const;

  virtual const std::string addToParameterDisplay() const {return "";}

};


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

  void  unaryIter(                           DensityOperatorLow&,  UnaryFunction) const;
  void binaryIter(const DensityOperatorLow&, DensityOperatorLow&, BinaryFunction) const;

  const std::string addToParameterDisplay() const {return " Fast Iteration.";}

  const blitzplusplus::SlicesData<2*RANK,blitzplusplus::vfmsi::LeftRight<RANK,blitzplusplus::vfmsi::Left> > slicesData_;

};



} // master



#define BASE_class boost::mpl::if_c<IS_FAST,master::BaseFast<RANK>,master::Base<RANK> >::type

/// An \link trajectory::Adaptive Adaptive\endlink trajectory class representing Master equation evolution from a \link quantumdata::DensityOperator density-operator\endlink initial condition
/**
 * Cf. \ref masterequation.
 *
 * \tparam RANK arity of the Hilbert space
 * \tparam V has the same function as the template parameter `V` in display_densityoperator::_, which class is used here for displaying the evolved density operator
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

  using Base::getTime; using Base::getAv;

  template<typename SYS>
  Master(DensityOperator& rho, const SYS& sys, const master::Pars& pt, bool negativity,
         const DensityOperatorLow& scaleAbs=DensityOperatorLow())
    : Base(rho,cpputils::sharedPointerize(sys),pt,scaleAbs), doDisplay_(getAv(),negativity)
  {}
  
  const DensityOperator& getRho() const {return rho_;}

private:
  using Base::rho_;

  std::ostream& display_v   (std::ostream& os, int precision) const {return doDisplay_.display   (getTime(),rho_,os,precision);}
  std::ostream& displayKey_v(std::ostream& os, size_t& i    ) const {return doDisplay_.displayKey(os,i);}
  const display_densityoperator::_<RANK,V> doDisplay_;

  const std::string trajectoryID_v() const {return "Master";}

};



} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MASTER_H_INCLUDED
