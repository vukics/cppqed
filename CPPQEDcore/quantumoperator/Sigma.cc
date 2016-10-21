// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#if !BOOST_PP_IS_ITERATING

#include "Sigma.h"

#include "TMP_Tools.h"

using quantumdata::Types;
using blitz::Range;

#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>

#define BOOST_PP_ITERATION_LIMITS (2,BOOST_PP_SUB(BLITZ_ARRAY_LARGEST_RANK,1))
#define BOOST_PP_FILENAME_1 "../quantumoperator/Sigma.cc"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

#else  // BOOST_PP_IS_ITERATING

#define ITER   BOOST_PP_ITERATION()
#define ITERM1 BOOST_PP_SUB(ITER,1)

#define PartialProject_print(z,nouse,unused) Range::all()

namespace quantumoperator {


template<>
const Types<ITERM1>::StateVectorLow
partialProject<ITER,true >(const Types<ITER>::StateVectorLow& psi, int n)
{
  return psi(n,BOOST_PP_ENUM(ITERM1,PartialProject_print,~));
}


template<>
const Types<ITERM1>::StateVectorLow
partialProject<ITER,false>(const Types<ITER>::StateVectorLow& psi, int n)
{
  return psi(BOOST_PP_ENUM(ITERM1,PartialProject_print,~),n);
}


} // quantumoperator

#undef PartialProject_print

#undef ITERM1
#undef ITER


#endif // BOOST_PP_IS_ITERATING
