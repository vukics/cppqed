// -*- C++ -*-
#ifndef   SIGMA_TYPE_OPERATOR_IMPL_INCLUDED
#define   SIGMA_TYPE_OPERATOR_IMPL_INCLUDED

#include "Sigma.h"

namespace quantumoperator {


template<int L, int R, typename OTHER, bool IS_HEAD>
class DirectProduct
{
public:
  static const int N_RANK=OTHER::N_RANK+1;

  typedef typename quantumdata::Types<N_RANK>::StateVectorLow StateVectorLow;

  DirectProduct(const OTHER& other) : other_(other) {}

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

#endif // SIGMA_TYPE_OPERATOR_IMPL_INCLUDED
