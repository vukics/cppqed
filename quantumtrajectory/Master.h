// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_MASTER_H_INCLUDED
#define QUANTUMTRAJECTORY_MASTER_H_INCLUDED

#include "MasterFwd.h"

#include "ExactFwd.h"
#include "HamiltonianFwd.h"
#include "LiouvilleanFwd.h"

#include "DO_Display.h"
#include "Types.h"

#include "Exception.h"
#include "SmartPtr.h"
#include "VectorFromMatrixSliceIterator.h"

#include <boost/function.hpp>



namespace quantumtrajectory {


namespace master {


typedef trajectory::ParsTrajectory Pars;


struct NonUnitaryIP  : cpputils::Exception {};
struct NoLiouvillean : cpputils::Exception {};


template<int RANK>
class Base : public trajectory::Trajectory<typename quantumdata::Types<RANK>::DensityOperatorLow>
{
public:
  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Exact        <RANK> Exact        ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian  ;
  typedef structure::Liouvillean  <RANK> Liouvillean  ;
  typedef structure::Averaged     <RANK> Averaged   ;

  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;
  typedef typename quantumdata::Types<RANK>::    StateVectorLow     StateVectorLow;

  typedef trajectory::Trajectory<DensityOperatorLow> TrajectoryBase;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;
  
  typedef structure::QuantumSystemWrapper<RANK> QuantumSystemWrapper;

  using TrajectoryBase::getEvolved; using TrajectoryBase::getDtDid; using TrajectoryBase::getTime; using TrajectoryBase::getOstream;

  Base(DensityOperator&, typename QuantumSystem::Ptr, const Pars&, const DensityOperatorLow& =DensityOperatorLow());

  void derivs(double, const DensityOperatorLow&, DensityOperatorLow&) const;

  void step             (double) const;

  void displayParameters(      ) const;

protected:
  typedef boost::function<void(                       StateVectorLow&)>  UnaryFunction;
  typedef boost::function<void(const StateVectorLow&, StateVectorLow&)> BinaryFunction;

  DensityOperator& rho_;

  const typename Averaged::Ptr getAv() const {return qs_.getAv();}

private:
  virtual void  unaryIter(                           DensityOperatorLow&,  UnaryFunction) const;
  virtual void binaryIter(const DensityOperatorLow&, DensityOperatorLow&, BinaryFunction) const;

  virtual const std::string addToParameterDisplay() const {return "";}

  // ****************
  // *** Data members

  mutable double tIntPic0_; // The time instant of the beginning of the current time step.

  const QuantumSystemWrapper qs_;

};


template<int RANK>
class BaseFast : public Base<RANK>
{
public:
  typedef typename Base<RANK>::QuantumSystem QuantumSystem;

  typedef typename Base<RANK>::DensityOperator DensityOperator;

  typedef typename Base<RANK>::DensityOperatorLow DensityOperatorLow;

  BaseFast(DensityOperator& rho, typename QuantumSystem::Ptr sys, const Pars& p, const DensityOperatorLow& scaleAbs=DensityOperatorLow())
    : Base<RANK>(rho,sys,p,scaleAbs), slicesData_(rho()) {}

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


template<int RANK, typename V, bool IS_FAST>
class Master : public BASE_class
{
public:
  typedef typename BASE_class Base;

#undef  BASE_class

  typedef typename Base::QuantumSystem QuantumSystem;

  typedef typename Base::DensityOperatorLow DensityOperatorLow; 

  typedef typename Base::DensityOperator DensityOperator;

  using Base::getOstream; using Base::getPrecision; using Base::getTime; using Base::getAv;

  template<typename SYS>
  Master(DensityOperator& rho, const SYS& sys, const master::Pars& pt, bool negativity,
	 const DensityOperatorLow& scaleAbs=DensityOperatorLow())
    : trajectory::TrajectoryBase(pt), Base(rho,cpputils::sharedPointerize(sys),pt,scaleAbs), doDisplay_(getAv(),pt,negativity)
  {}

private:
  using Base::rho_;

  void   displayMore   () const {doDisplay_.displayMore(getTime(),rho_,getOstream(),getPrecision());}
  size_t displayMoreKey() const {return doDisplay_.displayMoreKey(getOstream());}
  
  const details::DO_Display<RANK,V> doDisplay_;

};



} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MASTER_H_INCLUDED
