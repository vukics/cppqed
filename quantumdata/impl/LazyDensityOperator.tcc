// -*- C++ -*-
#ifndef   QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED
#define   QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED

#include "LazyDensityOperator.h"

#include "impl/LazyDensityOperatorSliceIterator.tcc"

namespace quantumdata {

template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>& matrix, F function)
{
  ldo::DiagonalIterator<RANK,V> begin(matrix.template begin<V>());
  T init(function(*begin));

  using namespace cpputils;
  return T(
           accumulate(++begin,
                      matrix.template end<V>(),
                      init,
                      function,
                      cpputils::plus<T>()
                      )
           );
}


template<int RANK> template<typename V>
const ldo::DiagonalIterator<RANK,V>
LazyDensityOperator<RANK>::begin() const
{
  return ldo::DiagonalIterator<RANK,V>(*this,mpl::false_());
}

template<int RANK> template<typename V>
const ldo::DiagonalIterator<RANK,V>
LazyDensityOperator<RANK>::end  () const
{
  return ldo::DiagonalIterator<RANK,V>(*this,mpl:: true_());
}


} //quantumdata

#endif // QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED
