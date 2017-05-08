// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZARRAYEXTENSIONS_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZARRAYEXTENSIONS_TCC_INCLUDED

#include "BlitzArrayExtensions.h"

#include "TMP_Tools.h"
#include "Conversions.h"

#include "BlitzTinyExtensions.tcc"


namespace blitzplusplus {

using blitz::shape;
using blitz::neverDeleteData;


template<typename T, int RANK>
const blitz::Array<T,1>
unaryArray(const blitz::Array<T,RANK>& array)
{
  if (!array.data()) return blitz::Array<T,1>();
#ifndef   NDEBUG
  if (!array.isStorageContiguous()) throw (NonContiguousStorageException());
#endif // NDEBUG
  return blitz::Array<T,1>(const_cast<T*>(array.data()),shape(array.size()),neverDeleteData);
}


template<typename T, int TWO_TIMES_RANK>
const blitz::Array<T,2>
binaryArray(const blitz::Array<T,TWO_TIMES_RANK>& array)
{
  if (!array.data()) return blitz::Array<T,2>();

#ifndef   NDEBUG
  if (!array.isStorageContiguous()) throw (NonContiguousStorageException());

  static const int RANK=tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>::value;

  {
    const blitz::TinyVector<int,TWO_TIMES_RANK>& ordering=array.ordering();
    bool correct=true;
    for (int i=0; i<RANK; i++)
      correct&=(ordering(i)==ordering(i+RANK)+RANK);
    if (!correct) throw BinaryArrayOrderingErrorException();
  }
#endif // NDEBUG

  int size=long2Int(product(halfCutTiny(array.shape())));
  // apparently, the product is a long
  return blitz::Array<T,2>(const_cast<T*>(array.data()),shape(size,size),neverDeleteData);
}


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_BLITZARRAYEXTENSIONS_TCC_INCLUDED
