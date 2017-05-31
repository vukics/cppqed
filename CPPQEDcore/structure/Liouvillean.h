// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_LIOUVILLEAN_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_LIOUVILLEAN_H_INCLUDED

#include "LiouvilleanFwd.h"

#include "LiouvilleanAveragedCommon.h"
#include "TimeDependence.h"

#include "LazyDensityOperator.h"
#include "StateVector.h"


namespace structure {

struct SuperoperatorNotImplementedException : cpputils::Exception
{
  explicit SuperoperatorNotImplementedException(size_t m) : m_(m) {}
  
  size_t m_;
};

/// The interface every system having Liouvillean time-evolution must present towards the trajectory drivers
/**
 * The time-evolution must be Markovian where the Lindblad form of the Master equation is the most general one possible:
 * \f[\dot\rho=\frac1{i\hbar}\comm{H}\rho+\sum_m\lp J_m\rho J_m^\dag-\frac12\comm{J_m^\dag J_m}{\rho}_+\rp\f]
 * The class represents the set of \f$J_m\f$ operators (Lindblads or quantum jump operators) of arbitrary number,
 * either for Master equation or Monte Carlo wave-function evolution (in the latter case it calculates the jump rates as well).
 * 
 * \tparamRANK
 * 
 * \note No matter how complex the quantum system, the framework always assigns a unique ordinal to each jump corresponding to every subsystem
 * 
 * \note It is always possible to forgo the explicit calculation of certain jump rates because the rate can be calculated also on the basis of the Liouvillean::actWithJ function by the 
 * \link quantumtrajectory::MCWF_Trajectory MCWF stepper\endlink. The fact that such a fallback is desired can be signalled by setting a negative value for the rate of the given jump 
 * (“special jump”). \see \ref specialjump
 * 
 * \see The design is very similar to Exact, but here the choice is not between TwoTime/OneTime dependence, but OneTime/NoTime.
 * 
 */
template<int RANK>
class Liouvillean : public quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> >
{
public:
  static const int N_RANK=RANK;

  typedef boost::shared_ptr<const Liouvillean> Ptr;
  
private:
  typedef quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> > Base;

public:
  typedef quantumdata::StateVector<RANK> StateVector;
  
  typedef typename Base::    StateVectorLow     StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::DArray1D Rates; ///< The 1D real array for storing the jump rates

  virtual ~Liouvillean() {}
  
  /// Returns the set of jump rates \f$\bra{\Psi}J_m^\dagger J_m\ket{\Psi},\f$ where the Lindblads are in general time-dependent
  /** Simply redirects to LiouvilleanAveragedCommonRanked::average, so that this function does not appear in the interface for implementers. */
  const Rates rates(double t, const StateVector& psi) const {return Base::average(t,psi);}

  /// Performs the quantum jump operation \f$\ket\Psi\rightarrow J_m(t)\ket\Psi\f$
  void actWithJ(double t,            ///<[in] \f$t\f$
                StateVectorLow& psi, ///<[in/out] \f$\ket\Psi\f$
                size_t m             ///<[out] \f$m\f$
               ) const {return actWithJ_v(t,psi,m);}

  /// Calculates \f$\Lcal\rho=J_m(t)\rho J_m(t)^\dagger\f$ and adds it to `drhodt`
  void actWithSuperoperator(double t,                      ///<[in] time
                            const DensityOperatorLow& rho, ///<[in] density operator
                            DensityOperatorLow& drhodt,    ///<[in/out] density operator
                            size_t m                       ///<[in] ordinal of jump operator
                           ) const {actWithSuperoperator_v(t,rho,drhodt,m);}

private:
  virtual void actWithJ_v(double, StateVectorLow&, size_t) const = 0;
  virtual void actWithSuperoperator_v(double, const DensityOperatorLow&, DensityOperatorLow&, size_t m) const {throw SuperoperatorNotImplementedException(m);}

};


/// Implements the general Liouvillean interface by dispatching the two possible \link time::DispatcherIsTimeDependent time-dependence levels\endlink
/**
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT describes whether the \f$J_m\f$s are time-dependent. `true`: OneTime – `false`: NoTime
 * 
 */
template<int RANK, bool IS_TIME_DEPENDENT>
class LiouvilleanTimeDependenceDispatched : public Liouvillean<RANK>
{
public:
  typedef typename Liouvillean<RANK>::StateVectorLow      StateVectorLow     ;
  typedef typename Liouvillean<RANK>::DensityOperatorLow  DensityOperatorLow ;
  typedef typename Liouvillean<RANK>::LazyDensityOperator LazyDensityOperator;
  typedef typename Liouvillean<RANK>::Rates               Rates              ;
  
  typedef typename time::DispatcherIsTimeDependent<IS_TIME_DEPENDENT>::type Time;

private:
  void         actWithJ_v(double t, StateVectorLow& psi, size_t lindbladNo) const final {actWithJ_v(Time(t),psi,lindbladNo);}   ///< Redirects the virtual inherited from Liouvillean<RANK>
  const Rates   average_v(double t, const LazyDensityOperator&  matrix) const final {return rates_v(Time(t),matrix);} ///< Redirects the virtual inherited from LiouvilleanAveragedCommonRanked

  void actWithSuperoperator_v(double t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t m) const final {actWithSuperoperator_v(Time(t),rho,drhodt,m);}
  
  virtual void        actWithJ_v(Time, StateVectorLow&, size_t   ) const = 0;
  virtual const Rates    rates_v(Time, const LazyDensityOperator&) const = 0;

  virtual void actWithSuperoperator_v(Time, const DensityOperatorLow&, DensityOperatorLow&, size_t m) const {throw SuperoperatorNotImplementedException(m);}

};


/** \page specialjump The notion of “special jumps” in the framework
 * 
 * Implementers of Liouvillean (typically, certain free elements) may decide to forgo the explicit calculation of the jump rate corresponding to any of the Lindblad operators,
 * because the action of the Lindblad on a state vector contains enough information about the jump rate as well. This is signalled by a negative jump rate towards the framework.
 * 
 * The motivation for this is usually that the explicit calculation of a given jump rate is cumbersome and error-prone in some cases.
 * 
 * The quantumtrajectory::MCWF_trajectory class, when encountering a negative jump rate, will reach for the Liouvillean::actWithJ function to calculate the corresponding rate.
 * 
 * Assume that \f$J_m\f$ is the set of Lindblads, and for \f$m=\text{at}\f$ the jump rate is found negative. In this case, \f$\ket{\Psi_\text{at}}=J_\text{at}\ket\Psi\f$ is
 * calculated and cached by the trajectory driver (cf. quantumtrajectory::MCWF_trajectory::calculateSpecialRates) and the given jump rate is taken as \f$\norm{\ket{\Psi_\text{at}}}^2\f$.
 * 
 * \note In certain cases, the use of a special jump can be less efficient than the explicit calculation of the rate.
 * 
 */


} // structure

#endif // CPPQEDCORE_STRUCTURE_LIOUVILLEAN_H_INCLUDED
