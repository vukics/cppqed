// -*- C++ -*-
#ifndef   _BLITZ_TO_FLENS_FWD_H
#define   _BLITZ_TO_FLENS_FWD_H

#include<flens/storage.h>


namespace blitz2flens {


using namespace flens;


template<typename>
struct DenseVectorMF;

template<typename, StorageOrder>
struct GeMatrixMF   ;


template<StorageOrder>
struct StorageTag {};


} // blitz2flens


#endif // _BLITZ_TO_FLENS_FWD_H
