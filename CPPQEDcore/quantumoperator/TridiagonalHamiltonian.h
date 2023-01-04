// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONALHAMILTONIAN_H_INCLUDED
#define CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONALHAMILTONIAN_H_INCLUDED

#include "Hamiltonian.h"

#include "Tridiagonal.tcc"

#include <list>


namespace quantumoperator {


namespace details {

using namespace structure;

/// TridiagonalHamiltonian base
template<int RANK, bool IS_TIME_DEPENDENT>
class TDH_Base 
  : public HamiltonianTimeDependenceDispatched<RANK,
                                               std::conditional_t<IS_TIME_DEPENDENT,
                                                                  tmptools::integral_c<TimeDependence::ONE>,
                                                                  tmptools::integral_c<TimeDependence::NO>
                                                                  >::value
                                               >
{
protected:
  typedef quantumoperator::Tridiagonal<RANK> Tridiagonal;
  typedef std::list<Tridiagonal>             Tridiagonals;

  template<typename... ArgumentPack>
  TDH_Base(ArgumentPack&&... a) : hOverIs_(a...) {}

  const Tridiagonals& getH_OverIs() const {return hOverIs_;}
        Tridiagonals& getH_OverIs()       {return const_cast<Tridiagonals&>(static_cast<const TDH_Base*>(this)->getH_OverIs());}

private:
  /// \name Virtual inherited from HamiltonianTimeDependenceDispatched
  void addContribution_v(std::conditional_t<IS_TIME_DEPENDENT,OneTime,NoTime> t, const StateVectorLow<RANK>& psi, StateVectorLow<RANK>& dpsidt) const override
  {
    for (auto& h : hOverIs_)
      if constexpr (IS_TIME_DEPENDENT) {apply<RANK>(psi,dpsidt,h.propagate(double(t)));}
      else {apply<RANK>(psi,dpsidt,h);}
  }

  mutable Tridiagonals hOverIs_;

};


} // details



/// Implements the action of a Hamiltonian whose matrix consists of a sum of \link Tridiagonal tridiagonal matrices\endlink
/**
 * \f[H_\text{tridiagonals}(t)=H_0(t)+H_1(t)+H_2(t)+\dots\f] with the \f$H_i(t)\f$s being all described by Tridiagonal `<RANK>` objects.
 * 
 * Such a class can be constructed with either a list of Tridiagonal `<RANK>` objects, or only one such object when the above sum consists of only one term.
 * 
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT governs time-dependence & the composition of the class @ compile time
 * 
 * Implements Hamiltonian `<RANK,ONE>` when `IS_TIME_DEPENDENT=true` **OR** Hamiltonian `<RANK,NO>`  when `IS_TIME_DEPENDENT=false`
 * 
 * \note The present architecture of Tridiagonal does not allow to cover the case structure::TWO.
 *  
 */
template<int RANK, bool IS_TIME_DEPENDENT>
class TridiagonalHamiltonian : public details::TDH_Base<RANK,IS_TIME_DEPENDENT>
{
private:
  typedef details::TDH_Base<RANK,IS_TIME_DEPENDENT> Base;

public:
  typedef typename Base::Tridiagonal  Tridiagonal ;
  typedef typename Base::Tridiagonals Tridiagonals;

  template<typename... ArgumentPack>
  TridiagonalHamiltonian(ArgumentPack&&... a) : Base(a...) {} ///< Generic constructor

  TridiagonalHamiltonian(Tridiagonal&& hOverI ) : TridiagonalHamiltonian(1,hOverI) {} ///< \overload

};

#undef BASE_class

} // quantumoperator


#endif // CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONALHAMILTONIAN_H_INCLUDED
