// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_TCC_INCLUDED

#include "LazyDensityOperator.h"

#include "LazyDensityOperatorSliceIterator.tcc"

#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/transformed.hpp>


namespace quantumdata {

template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>& matrix, F function)
{
  auto begin(matrix.template begin<V>());
  T init(function(*begin)); // we take a way around default constructing T, so that the implicit interface is less stringently defined

  using namespace boost;
  return accumulate(make_iterator_range(++begin,matrix.template end<V>()) | adaptors::transformed(function) , init );
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
const DArray<1> deflate(const LazyDensityOperator<RANK>& matrix, bool offDiagonals)
{
  using mathutils::sqr;
  
  const size_t dim=matrix.getTotalDimension();
  
  DArray<1> res(offDiagonals ? sqr(dim) : dim);
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(matrix.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  for (Iterator i(etalon); idx<dim; ++i)
    res(idx++)=matrix(*i);
  
  if (offDiagonals)
    for (Iterator i(etalon); idx<sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
        dcomp matrixElement(matrix(*i)(*j));
        res(idx++)=real(matrixElement);
        res(idx++)=imag(matrixElement);
      }

  return res;

}


} //quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_TCC_INCLUDED
