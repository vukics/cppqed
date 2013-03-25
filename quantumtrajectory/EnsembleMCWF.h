// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
#define QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED

#include "EnsembleMCWFFwd.h"

#include "DO_Display.h"
#include "MCWF_Trajectory.h"
#include "DensityOperator.h"

#include "SmartPtr.h"


namespace quantumtrajectory {


////////////////////////////////
//
// Ensemble of MCWF trajectories
//
////////////////////////////////


namespace ensemblemcwf {


#define STATE_VECTORS(r) boost::ptr_vector<quantumdata::StateVector<r> >

#define BASE_class trajectory::Ensemble< quantumdata::DensityOperator<RANK>&, const quantumdata::StateVector<RANK>& >

template<int RANK>
class Base
  : private boost::base_from_member<STATE_VECTORS(RANK) >,
    public BASE_class
{
public:
  typedef STATE_VECTORS(RANK) StateVectors;

#undef  STATE_VECTORS

  typedef boost::base_from_member<StateVectors> StateVectorsBase;

  typedef BASE_class Ensemble;

#undef  BASE_class

  typedef MCWF_Trajectory<RANK> Adaptive;

  typedef typename Adaptive::StateVector    StateVector   ;
  typedef typename Adaptive::StateVectorLow StateVectorLow; 

  typedef typename structure::QuantumSystem<RANK>::Ptr QuantumSystemPtr;

  typedef typename Ensemble::Impl Trajectories;


  Base(
       const StateVector&,
       QuantumSystemPtr,
       const ParsMCWF&,
       const StateVectorLow& =StateVectorLow()
       );

protected:
  const QuantumSystemPtr getQS() const {return qs_;}

private:
  // static helpers to constructor
  static std::auto_ptr<StateVectors> stateVectors(const StateVector& psi, size_t nTraj);
  static std::auto_ptr<Trajectories> trajectories(StateVectors& psis, QuantumSystemPtr qs, const ParsMCWF& p, const StateVectorLow& scaleAbs);

  
  const typename Ensemble::TBA_Type getInitializedTBA_v() const {rho_()=0; return rho_;}

  mutable quantumdata::DensityOperator<RANK> rho_;

  const QuantumSystemPtr qs_;

};


} // ensemblemcwf


#define BASE_class ensemblemcwf::Base<RANK>


template<int RANK, typename V>
class EnsembleMCWF : public BASE_class
{
public:
  typedef BASE_class Base;

#undef  BASE_class

  typedef details::DO_Display<RANK,V> DO_Display;

  typedef typename Base::StateVectorLow StateVectorLow; 

  typedef typename Base::Ensemble Ensemble;

  typedef typename Base      ::    StateVector     StateVector;
  typedef typename DO_Display::DensityOperator DensityOperator;

  using Base::getQS; using Base::getTime; using Base::toBeAveraged;

  template<typename SYS>
  EnsembleMCWF(
               const StateVector& psi,
               const SYS& sys,
               const ParsMCWF& p,
               bool negativity,
               const StateVectorLow& scaleAbs=StateVectorLow()
               )
    : Base(psi,cpputils::sharedPointerize(sys),p,scaleAbs), doDisplay_(structure::qsa<RANK>(getQS()),p,negativity) {}

private:
  std::ostream& display_v   (std::ostream& os, int precision) const {return doDisplay_.display   (getTime(),toBeAveraged(),os,precision);}
  std::ostream& displayKey_v(std::ostream& os, size_t& i    ) const {return doDisplay_.displayKey(os,i);}

  const DO_Display doDisplay_;

};


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_ENSEMBLEMCWF_H_INCLUDED
