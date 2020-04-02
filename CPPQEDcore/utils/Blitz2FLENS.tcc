// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZ2FLENS_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZ2FLENS_TCC_INCLUDED

#include "Blitz2FLENS.h"

#include "BlitzTinyExtensions.tcc"

#include <flens/flens.cxx>


namespace blitz2flens {


template<typename T, int RANK>
const DenseVectorOf<T> vector(blitz::Array<T,RANK>& array)
{
  typedef DenseVectorOf<T> DenseVector;
  
  typedef typename DenseVector::Engine ArrayView;

  assert(array.isStorageContiguous());
  // NEEDS_WORK more assertions are needed to ensure that the data will indeed be interpreted in the correct way by FLENS

  return DenseVector(ArrayView(array.size(),array.dataFirst()));
}


template<StorageOrder SO, typename T, int TWO_TIMES_RANK>
const GeMatrixOf<T,SO> matrix(blitz::Array<T,TWO_TIMES_RANK>& array)
{
  tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>();

  typedef GeMatrixOf<T,SO> GeMatrix;
  
  typedef typename GeMatrix::Engine FullStorageView;

  assert(!array.size() || array.isStorageContiguous());
  // NEEDS_WORK assert(SO==RowMajor && ) ordering!!!

  size_t dim=product(blitzplusplus::halfCutTiny(array.extent()));
  // halfCutTiny will throw if the dimensions are not the same

  return GeMatrix(FullStorageView(dim,dim,array.dataFirst(),dim));

}


template<StorageOrder SO, int TWO_TIMES_RANK>
const HeMatrixOf<SO> hermitianMatrix(CArray<TWO_TIMES_RANK>& array)
{
  tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>();

  typedef HeMatrixOf<SO> HeMatrix;

  typedef typename HeMatrixOf<SO>::Engine FullStorageView;

  assert(!array.size() || array.isStorageContiguous());
  // NEEDS_WORK assert(SO==RowMajor && ) ordering!!!

  size_t dim=product(blitzplusplus::halfCutTiny(array.extent()));
  // halfCutTiny will throw if the dimensions are not the same

  return HeMatrix(FullStorageView(dim,dim,array.dataFirst(),dim,0,0),flens::Upper);

}


template<int TWO_TIMES_RANK>
const CArray<tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>::value> ev(CArray<TWO_TIMES_RANK> m)
{
  using namespace flens;
  GeMatrixOf<dcomp,ColMajor> a(matrix<ColMajor>(m));

  CArray<tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>::value> res(blitzplusplus::halfCutTiny(m.extent())); ///< The eigenvalues as Blitz++ array

  DenseVector<ArrayView<dcomp> >
    v(vector(res)); ///< The eigenvalues as FLENS array

  lapack::ev(false,false,a,v,a,a);

  return res;
}


} // blitz2flens


#endif // CPPQEDCORE_UTILS_BLITZ2FLENS_TCC_INCLUDED
