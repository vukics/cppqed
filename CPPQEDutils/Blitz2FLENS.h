// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED

#include "BlitzArray.h"
#include "TMP_Tools.h"

#include "BlitzTinyExtensions.tcc"


#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif // USE_CXXLAPACK

#include <flens/flens.cxx>


namespace blitz2flens {


using flens::StorageOrder, flens::RowMajor, flens::ColMajor;


template<typename T>
using DenseVectorOf=flens::DenseVector<flens::ArrayView<T> >;

template<typename T, StorageOrder SO>
using GeMatrixOf=flens::GeMatrix<flens::FullStorageView<T,SO> >;

template<StorageOrder SO>
using HeMatrixOf=flens::HeMatrix<flens::FullStorageView<dcomp,SO> >;



template<typename T, int RANK>
auto vector(blitz::Array<T,RANK>& array)
{
  assert(array.isStorageContiguous());
  // NEEDS_WORK more assertions are needed to ensure that the data will indeed be interpreted in the correct way by FLENS

  return DenseVectorOf<T>(typename DenseVectorOf<T>::Engine(array.size(),array.dataFirst()));
}


template<StorageOrder SO, typename T, int TWO_TIMES_RANK>
auto matrix(blitz::Array<T,TWO_TIMES_RANK>& array)
{
  tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>();

  assert(!array.size() || array.isStorageContiguous());
  // NEEDS_WORK assert(SO==RowMajor && ) ordering!!!

  size_t dim=product(blitzplusplus::halfCutTiny(array.extent()));
  // halfCutTiny will throw if the dimensions are not the same

  return GeMatrixOf<T,SO>(typename GeMatrixOf<T,SO>::Engine(dim,dim,array.dataFirst(),dim));

}


template<StorageOrder SO, int TWO_TIMES_RANK>
auto hermitianMatrix(CArray<TWO_TIMES_RANK>& array)
{
  tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>();

  assert(!array.size() || array.isStorageContiguous());
  // NEEDS_WORK assert(SO==RowMajor && ) ordering!!!

  size_t dim=product(blitzplusplus::halfCutTiny(array.extent()));
  // halfCutTiny will throw if the dimensions are not the same

  return HeMatrixOf<SO>(typename HeMatrixOf<SO>::Engine(dim,dim,array.dataFirst(),dim,0,0),flens::Upper);

}


template<int TWO_TIMES_RANK>
auto ev(CArray<TWO_TIMES_RANK> m)
{
  using namespace flens;
  
  GeMatrixOf<dcomp,ColMajor> a(matrix<ColMajor>(m));

  CArray<tmptools::AssertEvenAndDivideBy2_v<TWO_TIMES_RANK>> res(blitzplusplus::halfCutTiny(m.extent())); ///< The eigenvalues as Blitz++ array

  DenseVector<ArrayView<dcomp>> v(vector(res)); ///< The eigenvalues as FLENS array

  lapack::ev(false,false,a,v,a,a);

  return res;
}


} // blitz2flens


#endif // CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED
