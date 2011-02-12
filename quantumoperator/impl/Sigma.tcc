// -*- C++ -*-
#ifndef   SIGMA_TYPE_OPERATOR_IMPL_INCLUDED
#define   SIGMA_TYPE_OPERATOR_IMPL_INCLUDED


namespace quantumoperator {


template<int L, int R, typename OTHER, bool IS_HEAD>
class directProduct
{
public:
  directProduct(const OTHER& other) : other_(other) {}

  template<int RANK>
  void
  apply(const typename quantumdata::Types<RANK>::StateVectorLow& psi, typename quantumdata::Types<RANK>::StateVectorLow& dpsidt) const
  {
    typename quantumdata::Types<RANK-1>::StateVectorLow dpsidtProjected(partialProject<RANK,IS_HEAD>(dpsidt,L));
    other_.template apply<RANK-1>(partialProject<RANK,IS_HEAD>(psi,R),dpsidtProjected);
  }

private:
  const OTHER& other_;

};


template<int L, int R, typename OTHER>
const directProduct<L,R,OTHER,true >
operator*(const Sigma<L,R>&, const OTHER& other)
{
  return directProduct<L,R,OTHER,true >(other);
}


template<int L, int R, typename OTHER>
const directProduct<L,R,OTHER,false>
operator*(const OTHER& other, const Sigma<L,R>&)
{
  return directProduct<L,R,OTHER,false>(other);
}


template<int L1, int R1, int L2, int R2>
const directProduct<L1,R1,Sigma<L2,R2>,true>
operator*(const Sigma<L1,R1>&, const Sigma<L2,R2>& sigma)
{
  return directProduct<L1,R2,Sigma<L2,R2>,true>(sigma);
}

} // quantumoperator

#endif // SIGMA_TYPE_OPERATOR_IMPL_INCLUDED
