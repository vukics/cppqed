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


template<int RANK>
const TTD_DARRAY(1) deflate(const LazyDensityOperator<RANK>& matrix, bool offDiagonals)
{
  using mathutils::sqr;
  typedef typename LazyDensityOperator<RANK>::Dimensions Dimensions;
  
  const size_t dim=matrix.getTotalDimension();
  
  TTD_DARRAY(1) res(offDiagonals ? sqr(dim) : dim);
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(Dimensions(size_t(0)),matrix.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  for (Iterator i(etalon); idx<dim; ++i)
    res(idx++)=matrix(dispatchLDO_index(*i));
  
  if (offDiagonals)
    for (Iterator i(etalon); idx<sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
        dcomp matrixElement(matrix(dispatchLDO_index(*i),dispatchLDO_index(*j)));
        res(idx++)=real(matrixElement);
        res(idx++)=imag(matrixElement);
      }

  return res;

}


} //quantumdata

#endif // QUANTUMDATA_IMPL_LAZYDENSITYOPERATOR_TCC_INCLUDED
