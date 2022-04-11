// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_QUANTUMOPERATOR_SIGMA_H_INCLUDED
#define   CPPQEDCORE_QUANTUMOPERATOR_SIGMA_H_INCLUDED

#include "Types.h"


namespace quantumoperator {



/// Stateless class implementing the unary quantumoperator \f$\ket{L}\bra{R}\f$.
/**
 * The template parameters are eventually converted to runtime data when the Sigma::apply function is called. At that point, dimensionality errors are detected in debug mode.
 * 
 * \note Sigma should perhaps inherit from DimensionsBookkeeper as this would allow for an earlier detection of dimensionality errors.
 */
template<int L, int R>
class Sigma
{
public:
  static const int N_RANK=1;

  typedef quantumdata::StateVectorLow<1> StateVectorLow;

  /// Application of the operator on a state vector can be implemented trivially
  void
  apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
  {
    dpsidt(L)+=psi(R);
  }  

  /// Hermitian conjugation also trivial
  const Sigma<R,L> dagger() const {return Sigma<R,L>();}

};


/// Helper for DirectProduct::apply
/**
 * Calculates the state-vector slice \f[\ket{\Psi^{\avr{1,2,3,…,\text{rank-2},\text{rank-1}}}(\iota_0=n)}\in\bigotimes_{i=1,2,3,…,\text{rank-2},\text{rank-1}}\HSpace_i\f]
 * if `IS_HEAD=true` and \f[\ket{\Psi^{\avr{0,1,2,…,\text{rank-3},\text{rank-2}}}(\iota_\text{rank-1}=n)}\in\bigotimes_{i=0,1,2,…,\text{rank-3},\text{rank-2}}\HSpace_i\f]
 * if `IS_HEAD=false`. The code is automatically generated for all template-parameter combinations (`RANK` up to #BLITZ_ARRAY_LARGEST_RANK) via preprocessor metaprogramming.
 * Cf. quantumoperator/Sigma.cc
 * ~~~{.sh}
 * g++ -P -E -Iutils/ -Iquantumoperator/ -Iquantumdata/ quantumoperator/Sigma.cc | tail -n128
 * ~~~
 * 
 * \note It’s better to convert `n` into a runtime variable because then we can use complete specializations of this function. Eventually it has to be converted anyway into an index of `psi`.
 */
template<int RANK, bool IS_HEAD>
const quantumdata::StateVectorLow<RANK-1>
partialProject(const quantumdata::StateVectorLow<RANK>& psi, int n);


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

  typedef quantumdata::StateVectorLow<N_RANK> StateVectorLow;

  DirectProduct(const OTHER& other) : other_(other) {} ///< The object has to store a reference to the `OTHER` object to enable DirectProduct::apply as an operator
  
  /// Applies the clause on a state vector `psi`, adding the result to `dpsidt`, analogously to Tridiagonal::apply and structure::Hamiltonian::addContribution
  /**
   * It is implemented as taking the partial projection (quantumoperator::partialProject) of the state vectors according to `L` and `R`
   * (which at this point are converted into runtime data), and calling the `apply` function of the `OTHER` type.
   */
  void
  apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
  {
    auto dpsidtProjected(partialProject<N_RANK,IS_HEAD>(dpsidt,L));
    other_.apply(partialProject<N_RANK,IS_HEAD>(psi,R),dpsidtProjected);
  }

private:
  const OTHER& other_;

};


/// Direct-product operator for Sigma \related Sigma
/**
 * Together with its overloads, this operator aptly demonstrates the use of DirectProduct and the meaning of its template parameters.
 * 
 * Mixing (by composition from the right) version.
 * 
 * \tparam OTHER the other operator type to mix with
 * 
 */
template<int L, int R, typename OTHER>
const DirectProduct<L,R,OTHER,true >
operator*(const Sigma<L,R>&, const OTHER& other)
{
  return DirectProduct<L,R,OTHER,true >(other);
}


/// Direct-product operator for Sigma \related Sigma
/**
 * Together with its overloads, this operator aptly demonstrates the use of DirectProduct and the meaning of its template parameters.
 * 
 * Mixing (by composition from the left) version.
 * 
 * \tparam OTHER the other operator type to mix with
 * 
 */
template<int L, int R, typename OTHER>
const DirectProduct<L,R,OTHER,false>
operator*(const OTHER& other, const Sigma<L,R>&)
{
  return DirectProduct<L,R,OTHER,false>(other);
}


/// Direct-product operator for Sigma \related Sigma Homogeneous (non-mixing) version.
template<int L1, int R1, int L2, int R2>
const DirectProduct<L1,R1,Sigma<L2,R2>,true>
operator*(const Sigma<L1,R1>&, const Sigma<L2,R2>& sigma)
{
  return DirectProduct<L1,R2,Sigma<L2,R2>,true>(sigma);
}


} // quantumoperator

#endif // CPPQEDCORE_QUANTUMOPERATOR_SIGMA_H_INCLUDED
