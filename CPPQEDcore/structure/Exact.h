// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_EXACT_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_EXACT_H_INCLUDED

#include "TimeDependence.h"
#include "Types.h"

namespace structure {


/// The interface every system that needs transformation between two quantum mechanical pictures must present towards the trajectory drivers
/**
 * Experience shows that even when a system uses interaction picture (which is automatically the case if any of its subsystems does) – that is, part of its dynamics is solved exactly – 
 * it may still want to calculate the jump operators and quantum averages in the normal picture. (cf. Cases 1 & 3 \link TimeDependence above\endlink.)
 * This is useful e.g. to reuse the code written for the non-interaction-picture case.
 * 
 * In this case, the framework has to be provided with some means to transform between the two pictures.
 * This is fulfilled by this class, from which classes describing such systems have to inherit.
 * 
 * E.g. if quantumtrajectory::MCWF_Trajectory sees that the simulated system inherits from Exact, then it will make the coherent part of the evolution in interaction picture,
 * whereupon it transforms back to normal picture, so that all the rest (jump rates, eventual jumps, quantum averages) can be calculated in this latter picture.
 * This makes that the two pictures coincide before each timestep. (Cf. also the stages described @ quantumtrajectory::MCWF_trajectory.)
 * 
 * \tparamRANK
 * \tparam IS_TWO_TIME default `true`, the most general case
 * 
 * \see The design is very similar to that of Hamiltonian & HamiltonianTimeDependenceDispatched.
 * 
 */
template<int RANK>
class Exact
{
public:
  virtual ~Exact() {}

  /// Describes the operation which transforms from interaction picture to the normal picture: \f$\ket{\Psi(t)}\rightarrow U(t,t_0)\ket{\Psi}\f$
  void actWithU(double t, ///<[in] \f$t\f$
                StateVectorLow<RANK>& psi, ///<[in/out] \f$\ket{\Psi}\f$
                double t0 ///<[in] \f$t_0\f$
               ) const {return actWithU_v(t,psi,t0);} 

private:
  virtual void actWithU_v(double, StateVectorLow<RANK>&, double) const = 0;

};


template <int RANK>
using ExactPtr=std::shared_ptr<const Exact<RANK>>;


/// Implements the general Exact interface by dispatching the two possible \link time::DispatcherIsTwoTime time-dependence levels\endlink
/**
 * \tparamRANK
 * \tparam IS_TWO_TIME `true`: TwoTime – `false`: OneTime
 * 
 */
template<int RANK, bool IS_TWO_TIME>
class ExactTimeDependenceDispatched : public Exact<RANK>
{
public:
  typedef time::DispatcherIsTwoTime_t<IS_TWO_TIME> Time;
  
private:
  void actWithU_v(double t, StateVectorLow<RANK>& psi, double t0) const final {actWithU_v(Time(t,t0),psi);}

  virtual void actWithU_v(Time, StateVectorLow<RANK>& psi) const = 0;

};


} // structure

#endif // CPPQEDCORE_STRUCTURE_EXACT_H_INCLUDED
