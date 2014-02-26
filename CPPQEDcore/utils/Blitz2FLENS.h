// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED

#include "BlitzArray.h"

#include <flens/flens.h>

#include <blitz/array.h>

#include <boost/mpl/identity.hpp>
namespace mpl=boost::mpl;


namespace blitz2flens {


using namespace flens;


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


template<StorageOrder SO, typename T, int TWO_TIMES_RANK>
const typename GeMatrixMF<T,SO>::type matrix(const blitz::Array<T,TWO_TIMES_RANK>&);



template<StorageOrder SO, int TWO_TIMES_RANK>
const typename HeMatrixMF<SO>::type hermitianMatrix(const CArray<TWO_TIMES_RANK>&);


} // blitz2flens

#endif // CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED
