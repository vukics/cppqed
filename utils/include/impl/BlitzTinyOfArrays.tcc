// -*- C++ -*-
#ifndef   _BLITZ_TINY_OF_ARRAYS_IMPL_H
#define   _BLITZ_TINY_OF_ARRAYS_IMPL_H

#include "BlitzTinyOfArrays.h"

#include "Range.h"

#include<boost/lambda/lambda.hpp>



// The following are implemented as runtime loops because they would be too slow at compile time

namespace blitzplusplus {


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(TOA_ShallowCopy, const T_numtype& initValue)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(initValue);
}


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(TOA_DeepCopy   , const T_numtype& initValue)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(initValue.copy());
}


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(TOA_ShallowCopy, const TinyOfArrays& x)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(x(i));
}


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(TOA_DeepCopy   , const TinyOfArrays& x)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(x(i).copy());
}



template<typename T, int RANK, int LENGTH> 
const TinyOfArrays<T,RANK,LENGTH>
negate(const TinyOfArrays<T,RANK,LENGTH>& arrays)
{
  TinyOfArrays<T,RANK,LENGTH> res(TOA_DeepCopy(),arrays);
  boost::for_each(res,boost::lambda::_1*=-1);
  return res;
}


} // blitzplusplus

#endif // _BLITZ_TINY_OF_ARRAYS_IMPL_H
