// -*- C++ -*-
#if !BOOST_PP_IS_ITERATING

#ifndef   UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED

#include "ComplexArrayExtensions.h"

#include "impl/BlitzArrayExtensions.tcc"

#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/arithmetic/div.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>


namespace blitzplusplus {


namespace details {

#define BOOST_PP_ITERATION_LIMITS (1,BOOST_PP_DIV(BLITZ_ARRAY_LARGEST_RANK,2))
#define BOOST_PP_FILENAME_1 "impl/ComplexArrayExtensions.tcc"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

} // details

template<int TWO_TIMES_RANK>
inline
void
hermitianConjugateSelf(CArray<TWO_TIMES_RANK>& array) 
{
  BOOST_STATIC_ASSERT( tmptools::IsEvenAssert<TWO_TIMES_RANK>::value );
  array=conj(array); 
  details::helpHermitianConjugate(array);
}


template<int TWO_TIMES_RANK>
inline
const CArray<TWO_TIMES_RANK>
hermitianConjugate(const CArray<TWO_TIMES_RANK>& array)
{
  BOOST_STATIC_ASSERT( tmptools::IsEvenAssert<TWO_TIMES_RANK>::value );
  CArray<TWO_TIMES_RANK> res(conj(array));
  details::helpHermitianConjugate(res);
  return res;
}


namespace dodirect {

using namespace linalg;

template<bool> void doDirect(CMatrix&, const CVector&, const CVector&);

const bool multiplication=true;
const bool addition=false;


} // dodirect



template<bool IS_MULTIPLICATION, int RANK1, int RANK2>
const CArray<RANK1+RANK2>
doDirect(const CArray<RANK1>& array1, const CArray<RANK2>& array2)
{
  using namespace linalg;
  CArray<RANK1+RANK2> res(concatenateTinies(array1.shape(),array2.shape()));
  if (res.data()) {
    CVector r1array1(unaryArray(array1));
    CVector r1array2(unaryArray(array2));
    CMatrix r2res(res.data(),shape(array1.size(),array2.size()),neverDeleteData);
    
    dodirect::doDirect<IS_MULTIPLICATION>(r2res,r1array1,r1array2);
  }
  return res;

}



} // blitzplusplus

#endif // UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED


#else  // BOOST_PP_IS_ITERATING


#define ITER BOOST_PP_ITERATION()
#define ITERT2 BOOST_PP_MUL(2,ITER)

#define HCH_print(z,m,unused) m ,

inline void 
helpHermitianConjugate(CArray<ITERT2>& array)
{
  array.transposeSelf(BOOST_PP_REPEAT_FROM_TO(ITER,ITERT2,HCH_print,~) BOOST_PP_ENUM_PARAMS(ITER,));
}

#undef HCH_print

#undef ITERT2
#undef ITER


#endif // BOOST_PP_IS_ITERATING
