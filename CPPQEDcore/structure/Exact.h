// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_EXACT_H_INCLUDED
#define STRUCTURE_EXACT_H_INCLUDED

#include "ExactFwd.h"

#include "Time.h"
#include "Types.h"

#include <boost/shared_ptr.hpp>

namespace structure {


/// The template-parameter-independent base of Exact
class ExactCommon
{
public:
  typedef boost::shared_ptr<const ExactCommon> Ptr;

  virtual ~ExactCommon() {}

  bool isUnitary() const {return isUnitary_v();} ///< Describes whether the interaction picture is unitary

private:
  virtual bool isUnitary_v() const = 0;

};


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
class Exact : public ExactCommon, private quantumdata::Types<RANK> 
{
public:
  typedef boost::shared_ptr<const Exact> Ptr;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  virtual ~Exact() {}

  /// Describes the operation which transforms from interaction picture to the normal picture: \f$\ket{\Psi(t)}\rightarrow U(t,t_0)\ket{\Psi}\f$
  void actWithU(double t, ///<[in] \f$t\f$
                StateVectorLow& psi, ///<[in/out] \f$\ket{\Psi}\f$
                double tIntPic0 ///<[in] \f$t_0\f$
               ) const {return actWithU_v(t,psi,tIntPic0);} 

private:
  virtual void actWithU_v(double, StateVectorLow&, double) const = 0;

};


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
  typedef typename Exact<RANK>::StateVectorLow StateVectorLow;
  
  typedef typename time::DispatcherIsTwoTime<IS_TWO_TIME>::type Time;
  
private:
  void actWithU_v(double t, StateVectorLow& psi, double tIntPic0) const {actWithU_v(Time(t,tIntPic0),psi);}

  virtual void actWithU_v(Time, StateVectorLow& psi) const = 0;

};


} // structure

#endif // STRUCTURE_EXACT_H_INCLUDED
