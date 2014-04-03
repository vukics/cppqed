// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_HAMILTONIAN_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_HAMILTONIAN_H_INCLUDED

#include "HamiltonianFwd.h"

#include "TimeDependence.h"
#include "Types.h"

#include <boost/shared_ptr.hpp>


namespace structure {


/// The interface every system having (possibly non-Hermitian) Hamiltonian time-evolution must present towards the trajectory drivers
/**
 * \tparamRANK
 *
 * \todo Hamiltonian should be able not only to add its contribution to a state vector, but simply return its contribution.
 * In Hamiltonian::addContribution, somehow signal to the function whether it has to add or replace, because replacement could be also useful in some contexts.
 *
 */
template<int RANK>
class Hamiltonian : private quantumdata::Types<RANK>
{
public:
  typedef boost::shared_ptr<const Hamiltonian> Ptr;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  /// Adds the Hamiltonian contribution \f$\frac{H(t)}i\ket\Psi\f$ of the given (sub)system to `dpsidt`
  /**
   * The assumption is that the time when the Schrödinger picture and interaction picture (if any) coincide is `t0`. There are two important points to note:
   * -# The contribution has to be *added* to `dpsidt` instead of `dpsidt` being *replaced*. 
   * This is because when the given system is embedded in a larger system, other (sub)systems may also contribute.
   * -# The function has to calculate the effect of \f$\frac{H(t)}i\f$ and not merely \f$H\f$, since it is the former which determines the derivative of the state vector.
   * 
   * This latter is so often missed, that we emphasize it again (although we know that it will still be missed from time to time):
   * \warning When implementing the Hamiltonian, not \f$H\f$ itself but \f$\frac Hi\f$ has to supplied!
   * 
   * \todo 3rd argument should be passed by rvalue reference — move constructor for `blitz::Array` needed for this?
   *
   */
  void addContribution(double t, ///<[in] the time instant \f$t\f$ for #TWO_TIME dependence
                       const StateVectorLow& psi, ///<[in] the state vector \f$\ket\Psi\f$
                       StateVectorLow& dpsidt, ///<[in/out] the state vector to be contributed to by \f$\frac{H(t)}i\ket\Psi\f$
                       double t0 ///<[in] the time instant \f$t_0\f$ for #TWO_TIME dependence
                      ) const
                       {addContribution_v(t,psi,dpsidt,t0);}

  virtual ~Hamiltonian() {}

private:
  virtual void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const = 0; 

};

/// Implements the general Hamiltonian interface by dispatching the different \link time::Dispatcher time-dependence levels\endlink
/**
 * \tparamRANK
 * \tparam TD Degree of \link TimeDependence time dependence\endlink.
 * 
 */
template<int RANK, TimeDependence TD>
class HamiltonianTimeDependenceDispatched : public Hamiltonian<RANK>
{
public:
  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow;
  
  typedef typename time::Dispatcher<TD>::type Time; ///< The actual time-dependence level from the template parameter `TD`

private:
  /// The inherited virtual gets implemented by calling a newly defined virtual with only one parameter describing different \link time::Dispatcher time-dependence levels\endlink
  void addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const final
  {
    addContribution_v(Time(t,t0),psi,dpsidt);
  }

  /// The newly defined virtual with only one parameter describing different \link time::Dispatcher time-dependence levels\endlink
  virtual void addContribution_v(Time, const StateVectorLow&, StateVectorLow&) const = 0;
  
};


} // structure

#endif // CPPQEDCORE_STRUCTURE_HAMILTONIAN_H_INCLUDED
