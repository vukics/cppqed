// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines quantumoperator::DirectProduct & helpers}
#ifndef   CPPQEDCORE_QUANTUMOPERATOR_SIGMA_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMOPERATOR_SIGMA_TCC_INCLUDED

#include "Sigma.h"

namespace quantumoperator {


/// A direct-product clause, representing a germinal expression-template mechanism for direct products of Sigma`<L,R>` instances with `OTHER` types.
/**
 * \tparam OTHER The type of the other member of the clause. Models:
 * -# Another class Sigma`<L1,R1>`
 * -# A Tridiagonal type
 * -# Another DirectProduct type (recursivity)
 * \tparam IS_HEAD Signifies whether Sigma`<L,R>` stands at the head or at the tail of the direct product
 * 
 * The class is nicely described by the signature of the related `operator*` operators
 */
template<int L, int R, typename OTHER, bool IS_HEAD>
class DirectProduct
{
public:
  static const int N_RANK=OTHER::N_RANK+1; ///< Reports the rank of the class for recursive usage

  typedef typename quantumdata::Types<N_RANK>::StateVectorLow StateVectorLow;

  DirectProduct(const OTHER& other) : other_(other) {} ///< The object has to store a reference to the `OTHER` object to enable DirectProduct::apply as an operator
  
  /// Applies the clause on a state vector `psi`, adding the result to `dpsidt`, analogously to Tridiagonal::apply and structure::Hamiltonian::addContribution
  /**
   * It is implemented as taking the partial projection (quantumoperator::partialProject) of the state vectors according to `L` and `R`
   * (which at this point are converted into runtime data), and calling the `apply` function of the `OTHER` type.
   */
  void
  apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
  {
    typename quantumdata::Types<N_RANK-1>::StateVectorLow dpsidtProjected(partialProject<N_RANK,IS_HEAD>(dpsidt,L));
    other_.apply(partialProject<N_RANK,IS_HEAD>(psi,R),dpsidtProjected);
  }

private:
  const OTHER& other_;

};


template<int L, int R, typename OTHER>
const DirectProduct<L,R,OTHER,true >
operator*(const Sigma<L,R>&, const OTHER& other)
{
  return DirectProduct<L,R,OTHER,true >(other);
}


template<int L, int R, typename OTHER>
const DirectProduct<L,R,OTHER,false>
operator*(const OTHER& other, const Sigma<L,R>&)
{
  return DirectProduct<L,R,OTHER,false>(other);
}


template<int L1, int R1, int L2, int R2>
const DirectProduct<L1,R1,Sigma<L2,R2>,true>
operator*(const Sigma<L1,R1>&, const Sigma<L2,R2>& sigma)
{
  return DirectProduct<L1,R2,Sigma<L2,R2>,true>(sigma);
}

} // quantumoperator

#endif // CPPQEDCORE_QUANTUMOPERATOR_SIGMA_TCC_INCLUDED
