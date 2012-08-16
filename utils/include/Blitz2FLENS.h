// -*- C++ -*-
#ifndef   UTILS_INCLUDE_BLITZ2FLENS_H_INCLUDED
#define   UTILS_INCLUDE_BLITZ2FLENS_H_INCLUDED

#include "BlitzArray.h"
#include "ComplexExtensions.h"

#include<flens/flens.h>

#include<blitz/array.h>

#include<boost/mpl/identity.hpp>
namespace mpl=boost::mpl;


namespace blitz2flens {


using namespace flens;

template<StorageOrder>
struct StorageTag {};

typedef StorageTag<RowMajor> RowMajorTag;
typedef StorageTag<ColMajor> ColMajorTag;


template<typename T>
struct DenseVectorMF : mpl::identity<DenseVector<ArrayView<T> > > 
{
  typedef ArrayView<T> View;
};


template<typename T, StorageOrder SO>
struct GeMatrixMF : mpl::identity<GeMatrix<FullStorageView<T,SO> > >
{
  typedef FullStorageView<T,SO> View;
};

template<StorageOrder SO>
struct HeMatrixMF : mpl::identity<HeMatrix<FullStorageView<dcomp,SO> > >
{
  typedef FullStorageView<dcomp,SO> View;
};


template<typename T, int RANK>
const typename DenseVectorMF<T>::type vector(const blitz::Array<T,          RANK>&);


template<typename T, int TWO_TIMES_RANK, StorageOrder SO>
const typename GeMatrixMF<T,SO>::type matrix(const blitz::Array<T,TWO_TIMES_RANK>&, StorageTag<SO> =StorageTag<SO>());



template<int TWO_TIMES_RANK, StorageOrder SO>
const typename HeMatrixMF<SO>::type hermitianMatrix(const TTD_CARRAY(TWO_TIMES_RANK)&, StorageTag<SO> =StorageTag<SO>());


} // blitz2flens

#endif // UTILS_INCLUDE_BLITZ2FLENS_H_INCLUDED
