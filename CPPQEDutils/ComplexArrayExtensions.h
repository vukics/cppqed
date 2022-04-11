// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Helpers for complex `blitz::Array`s, e.g. Hermitian conjugation of multi-matrices}
#if !BOOST_PP_IS_ITERATING

#ifndef   CPPQEDCORE_UTILS_COMPLEXARRAYEXTENSIONS_H_INCLUDED
#define   CPPQEDCORE_UTILS_COMPLEXARRAYEXTENSIONS_H_INCLUDED

#include "BlitzArray.h"
#include "BlitzArrayExtensions.h"
#include "CMatrix.h"
#include "MathExtensions.h"

#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/arithmetic/div.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>


namespace blitzplusplus {

using cppqedutils::sqrAbs;

BZ_DECLARE_FUNCTION_RET(sqrAbs,double)


namespace details {

#define BOOST_PP_ITERATION_LIMITS (1,BOOST_PP_DIV(BLITZ_ARRAY_LARGEST_RANK,2))
#define BOOST_PP_FILENAME_1 "ComplexArrayExtensions.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

} // details


template<int TWO_TIMES_RANK>
std::enable_if_t<TWO_TIMES_RANK%2==0>
hermitianConjugateSelf(CArray<TWO_TIMES_RANK>& array) 
{
  array=conj(array); 
  details::helpHermitianConjugate(array);
}


template<int TWO_TIMES_RANK>
std::enable_if_t<TWO_TIMES_RANK%2==0,CArray<TWO_TIMES_RANK>>
hermitianConjugate(const CArray<TWO_TIMES_RANK>& array)
{
  CArray<TWO_TIMES_RANK> res(conj(array));
  details::helpHermitianConjugate(res);
  return res;
}


namespace dodirect {

using namespace linalg;

template<bool> void _(CMatrix&, const CVector&, const CVector&);

const bool multiplication=true;
const bool addition=false;


} // dodirect


/// Direct product/sum
/**
 * Returns the direct product (if `IS_MULTIPLICATION`) \f$A_{i,j}=A1_i*A2_j\f$, or direct sum (otherwise)
 * \f$A_{i,j}=A1_i+A2_j\f$ of `array1` and `array2`, with \f$i,\ j\f$ running through all the multi-indices.
 * The implementation is in terms of blitzplusplus::unaryArray views of the arguments.
 * 
 * \tparam IS_MULTIPLICATION governs whether direct product or sum is to be calculated
*/
template<bool IS_MULTIPLICATION, int RANK1, int RANK2>
const CArray<RANK1+RANK2>
doDirect(const CArray<RANK1>& array1, const CArray<RANK2>& array2)
{
  using namespace linalg;
  CArray<RANK1+RANK2> res(concatenateTinies(array1.shape(),array2.shape()));
  if (res.data()) {
    CVector r1array1(unaryArray(array1));
    CVector r1array2(unaryArray(array2));
    CMatrix r2res(res.data(),blitz::shape(array1.size(),array2.size()),blitz::neverDeleteData);
    
    dodirect::_<IS_MULTIPLICATION>(r2res,r1array1,r1array2);
  }
  return res;

}


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_COMPLEXARRAYEXTENSIONS_H_INCLUDED


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
