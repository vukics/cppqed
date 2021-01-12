// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED

#include "Structure.h"

#include "StreamDensityOperator.h"
#include "QuantumTrajectory.h"
#include "Types.h"

#include "VectorFromMatrixSliceIterator.h"


namespace quantumtrajectory {


/// Auxiliary tools to Master
namespace master {


typedef trajectory::ParsEvolved Pars;


/// Thrown if the system is not applicable in Master-equation evolution
/**
 * \see structure::ExactCommon::applicableInMaster, \ref masterequationlimitations
 */
struct SystemNotApplicable : std::runtime_error {SystemNotApplicable() : std::runtime_error("") {}};


} // master



/// An \link trajectory::Adaptive Adaptive\endlink trajectory class representing Master equation evolution from a \link quantumdata::DensityOperator density-operator\endlink initial condition
/**
 * \see \ref masterequation.

 * \note The ODE driver underlying this class needs to store several (typically 6–7, this depends on the chosen driver) density-operator instants.
 *
 * \tparam RANK arity of the Hilbert space
 * \tparam V has the same function as the template parameter `V` in stream_densityoperator::_, which class is used here for deriving quantum averages to stream from the evolved density operator
 * 
 */
template<int RANK, typename V=tmptools::V_Empty>
class Master : public QuantumTrajectory<RANK,
                                        trajectory::Adaptive<typename quantumdata::Types<RANK>::DensityOperatorLow,
                                                             trajectory::Trajectory<typename structure::AveragedCommon::Averages>>>
{
public:
  using Adaptive=trajectory::Adaptive<typename quantumdata::Types<RANK>::DensityOperatorLow, trajectory::Trajectory<typename structure::AveragedCommon::Averages>>;
  
  using QTraj=QuantumTrajectory<RANK,Adaptive>;
  
  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Exact        <RANK> Exact        ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian  ;
  typedef structure::Liouvillean  <RANK> Liouvillean  ;
  typedef structure::Averaged     <RANK> Averaged     ;

  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;
  typedef typename quantumdata::Types<RANK>::    StateVectorLow     StateVectorLow;

  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  /// The actual function calculating the time derivative for \link evolved::Evolved ODE evolution\endlink
  void derivs(double, const DensityOperatorLow&, DensityOperatorLow&) const;

  Master(DensityOperator&& rho, ///< the density operator to be evolved
         typename QuantumSystem::Ptr sys, ///< object representing the quantum system
         const master::Pars& pt, ///< parameters of the evolution
         bool negativity, ///< governs whether entanglement should be calculated, cf. stream_densityoperator::_, quantumdata::negPT
         const DensityOperatorLow& scaleAbs=DensityOperatorLow() ///< has the same role as `scaleAbs` in evolved::Maker::operator()
        );

private:
  typename QTraj::StreamReturnType stream_v(std::ostream& os, int precision) const override {return dos_.stream(this->getTime(),rho_,os,precision);}
  
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return dos_.streamKey(os,i);}

  void step_v(double) final;

  std::ostream& streamParameters_v(std::ostream&) const override;

  const std::string trajectoryID_v() const override {return "Master";}

  DensityOperator rho_;

  const stream_densityoperator::_<RANK,V> dos_;

};



} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MASTER_H_INCLUDED
