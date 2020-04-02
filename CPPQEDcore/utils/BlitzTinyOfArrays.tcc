// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZTINYOFARRAYS_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZTINYOFARRAYS_TCC_INCLUDED

#include "BlitzTinyOfArrays.h"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/lambda/lambda.hpp>



// The following are implemented as runtime loops because they would be too slow at compile time

namespace blitzplusplus {


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(ShallowCopy, const T_numtype& initValue)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(initValue);
}


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(DeepCopy   , const T_numtype& initValue)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(initValue.copy());
}


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(ShallowCopy, const TinyOfArrays& x)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(x(i));
}


template<typename T, int RANK, int LENGTH>
inline TinyOfArrays<T,RANK,LENGTH>::TinyOfArrays(DeepCopy   , const TinyOfArrays& x)
{
  for (int i=0; i<LENGTH; ++i)
    (*this)(i).reference(x(i).copy());
}



template<typename T, int RANK, int LENGTH> 
const TinyOfArrays<T,RANK,LENGTH>
negate(const TinyOfArrays<T,RANK,LENGTH>& arrays)
{
  TinyOfArrays<T,RANK,LENGTH> res(DeepCopy(),arrays);
  boost::for_each(res,boost::lambda::_1*=-1);
  return res;
}


} // blitzplusplus

#endif // CPPQEDCORE_UTILS_BLITZTINYOFARRAYS_TCC_INCLUDED
