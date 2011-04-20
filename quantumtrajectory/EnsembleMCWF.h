// -*- C++ -*-
#ifndef _ENSEMBLE_OF_MCWF_TRAJECTORIES_H
#define _ENSEMBLE_OF_MCWF_TRAJECTORIES_H

#include "EnsembleMCWFFwd.h"

#include "DO_Display.h"
#include "MCWF_Trajectory.h"
#include "DensityOperator.h"


namespace quantumtrajectory {


////////////////////////////////
//
// Ensemble of MCWF trajectories
//
////////////////////////////////


namespace ensemblemcwf {


namespace details {

template<int RANK>
struct StateVectorsMF {

  typedef quantumdata::StateVector<RANK> StateVector;
  typedef boost::ptr_list<StateVector> type;

};

} // details


#define BASE_class trajectory::EnsembleTrajectories< quantumdata::DensityOperator<RANK>&, const quantumdata::StateVector<RANK>& >

template<int RANK>
class Base
  : private boost::base_from_member<typename details::StateVectorsMF<RANK>::type>,
    public BASE_class
{
public:
  typedef boost::base_from_member<typename details::StateVectorsMF<RANK>::type> StateVectorsBase;

  typedef BASE_class EnsembleTrajectories;

#undef  BASE_class

  typedef MCWF_Trajectory<RANK> Trajectory;

  typedef typename Trajectory::StateVector    StateVector   ;

  typedef typename Trajectory::StateVectorLow StateVectorLow; 

  typedef structure::QuantumSystem<RANK> QuantumSystem;

  typedef typename details::StateVectorsMF<RANK>::type StateVectors;
  typedef typename EnsembleTrajectories::Impl          Trajectories;


  Base(
       const StateVector&,
       const QuantumSystem&,
       const ParsMCWF_Trajectory&,
       const StateVectorLow& =StateVectorLow()
       );

private:
  mutable quantumdata::DensityOperator<RANK> rho_;

  const typename EnsembleTrajectories::TBA_Type getInitializedTBA() const {rho_()=0; return rho_;}

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

  typedef typename Base::QuantumSystem QuantumSystem;

  typedef typename Base::StateVectorLow StateVectorLow; 

  typedef typename Base::EnsembleTrajectories EnsembleTrajectories;

  typedef typename Base      ::    StateVector     StateVector;
  typedef typename DO_Display::DensityOperator DensityOperator;

  using Base::getOstream;

  EnsembleMCWF(
	       const StateVector& psi,
	       const QuantumSystem& sys,
	       const ParsMCWF_Trajectory& p,
	       bool negativity,
	       const StateVectorLow& scaleAbs=StateVectorLow()
	       )
    : trajectory::TrajectoryBase(p), Base(psi,sys,p,scaleAbs), doDisplay_(sys,p,negativity) {}

private:
  void   displayMore   (int precision) const {doDisplay_.displayMore(EnsembleTrajectories::toBeAveraged(),getOstream(),precision);}
  size_t displayMoreKey(             ) const {return doDisplay_.displayMoreKey(getOstream());}

  const DO_Display doDisplay_;

};


} // quantumtrajectory


#include "impl/EnsembleMCWF.tcc"


#endif // _ENSEMBLE_OF_MCWF_TRAJECTORIES_H
