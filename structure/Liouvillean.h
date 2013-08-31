// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_LIOUVILLEAN_H_INCLUDED
#define STRUCTURE_LIOUVILLEAN_H_INCLUDED

#include "LiouvilleanFwd.h"

#include "LiouvilleanAveragedCommon.h"

#include "LazyDensityOperator.h"
#include "StateVector.h"


namespace structure {


/// The first partial specialization of the general template Liouvillean for the one-time dependence case (\link TimeDependence Cases 1 & 2\endlink)
template<int RANK>
class Liouvillean<RANK,true> : public quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> >
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

private:
  virtual void actWithJ_v(double, StateVectorLow&, size_t) const = 0;

};


/// The second partial specialization of the general template Liouvillean for the no-time dependence case (\link TimeDependence Cases 3 & 4\endlink)
template<int RANK>
class Liouvillean<RANK,false> : public Liouvillean<RANK,true>
{
public:
  typedef typename Liouvillean<RANK,true>::StateVectorLow      StateVectorLow     ;
  typedef typename Liouvillean<RANK,true>::LazyDensityOperator LazyDensityOperator;
  typedef typename Liouvillean<RANK,true>::Rates               Rates              ;

private:
  void        actWithJ_v(double, StateVectorLow& psi, size_t jumpNo) const {actWithJ_v(psi,jumpNo);}   ///< Redirects the virtual inherited from Liouvillean<RANK,true>
  const Rates  average_v(double, const LazyDensityOperator&  matrix) const {return average_v(matrix);} ///< Redirects the virtual inherited from LiouvilleanAveragedCommonRanked

  virtual void        actWithJ_v(StateVectorLow&, size_t   ) const = 0;
  virtual const Rates  average_v(const LazyDensityOperator&) const = 0;

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

#endif // STRUCTURE_LIOUVILLEAN_H_INCLUDED
