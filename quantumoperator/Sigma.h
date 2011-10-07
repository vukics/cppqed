// -*- C++ -*-
#ifndef   SIGMA_TYPE_OPERATOR_INCLUDED
#define   SIGMA_TYPE_OPERATOR_INCLUDED

#include "SigmaFwd.h"

#include "Types.h"


namespace quantumoperator {


template<int L, int R, typename OTHER>
const directProduct<L,R,OTHER,true >
operator*(const Sigma<L,R>&, const OTHER&);

template<int L, int R, typename OTHER>
const directProduct<L,R,OTHER,false>
operator*(const OTHER&, const Sigma<L,R>&);

template<int L1, int R1, int L2, int R2>
const directProduct<L1,R1,Sigma<L2,R2>,true>
operator*(const Sigma<L1,R1>&, const Sigma<L2,R2>&);


template<int L, int R>
class Sigma
{
public:
  static const int N_RANK=1;

  typedef quantumdata::Types<1>::StateVectorLow StateVectorLow;

  void
  apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
  {
    dpsidt(L)+=psi(R);
  }  

  const Sigma<R,L> dagger() const {return Sigma<R,L>();}

};


template<int RANK, bool IS_HEAD>
const typename quantumdata::Types<RANK-1>::StateVectorLow
partialProject(const typename quantumdata::Types<RANK>::StateVectorLow&, int n);
// it's better to convert n into a runtime variable because then we can use complete specializations of this function


} // quantumoperator


#include "impl/Sigma.tcc"


#endif // SIGMA_TYPE_OPERATOR_INCLUDED
