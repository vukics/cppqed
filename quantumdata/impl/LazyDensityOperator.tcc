// -*- C++ -*-
#ifndef   QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED
#define   QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED

#include "LazyDensityOperator.h"

#include "impl/LazyDensityOperatorSliceIterator.tcc"

namespace quantumdata {

template<int RANK, typename F, typename V, typename T>
const T
partialTrace(const LazyDensityOperator<RANK>& matrix, F function, V v, T)
{
  ldo::DiagonalIterator<RANK,V> begin(matrix.begin(v));
  T init(function(*begin));

  using namespace cpputils;
  return T(
           accumulate(++begin,
                      matrix.end(v),
                      init,
                      function,
                      cpputils::plus<T>()
                      )
           );
}


template<int RANK> template<typename V>
const ldo::DiagonalIterator<RANK,V>
LazyDensityOperator<RANK>::begin(V) const
{
  return ldo::DiagonalIterator<RANK,V>(*this,mpl::false_());
}

template<int RANK> template<typename V>
const ldo::DiagonalIterator<RANK,V>
LazyDensityOperator<RANK>::end  (V) const
{
  return ldo::DiagonalIterator<RANK,V>(*this,mpl:: true_());
}


} //quantumdata

#endif // QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED
