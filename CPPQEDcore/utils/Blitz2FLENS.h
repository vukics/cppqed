// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED

#include "BlitzArray.h"
#include "TMP_Tools.h"

#ifndef USE_CXXLAPACK
#define USE_CXXLAPACK
#endif // USE_CXXLAPACK

#include <flens/flens.h>

#include <blitz/array.h>

#include <boost/mpl/identity.hpp>
namespace mpl=boost::mpl;


namespace blitz2flens {


using flens::StorageOrder; using flens::RowMajor; using flens::ColMajor;


template<typename T>
using DenseVectorOf=flens::DenseVector<flens::ArrayView<T> >;

template<typename T, StorageOrder SO>
using GeMatrixOf=flens::GeMatrix<flens::FullStorageView<T,SO> >;

template<StorageOrder SO>
using HeMatrixOf=flens::HeMatrix<flens::FullStorageView<dcomp,SO> >;


template<typename T, int RANK>
const DenseVectorOf<T> vector(const blitz::Array<T,          RANK>&);


template<StorageOrder SO, typename T, int TWO_TIMES_RANK>
const GeMatrixOf<T,SO> matrix(const blitz::Array<T,TWO_TIMES_RANK>&);


template<StorageOrder SO, int TWO_TIMES_RANK>
const HeMatrixOf<SO> hermitianMatrix(const CArray<TWO_TIMES_RANK>&);


template<int TWO_TIMES_RANK>
const CArray<tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>::value> ev(CArray<TWO_TIMES_RANK>);


} // blitz2flens

#endif // CPPQEDCORE_UTILS_BLITZ2FLENS_H_INCLUDED
