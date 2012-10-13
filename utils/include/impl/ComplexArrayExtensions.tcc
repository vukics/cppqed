// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED

#include "ComplexArrayExtensions.h"

#include "impl/BlitzArrayExtensions.tcc"
#include "impl/BlitzTinyExtensions.tcc"

#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/arithmetic/div.hpp>


namespace blitzplusplus {


namespace details {

#define BOOST_PP_ITERATION_LIMITS (1,BOOST_PP_DIV(BLITZ_ARRAY_LARGEST_RANK,2))
#define BOOST_PP_FILENAME_1 "details/HCH_ImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

} // details

template<int TWO_TIMES_RANK>
inline
void
hermitianConjugateSelf(TTD_CARRAY(TWO_TIMES_RANK)& array) 
{
  BOOST_STATIC_ASSERT( tmptools::IsEvenAssert<TWO_TIMES_RANK>::value );
  array=conj(array); 
  details::helpHermitianConjugate(array);
}


template<int TWO_TIMES_RANK>
inline
const TTD_CARRAY(TWO_TIMES_RANK)
hermitianConjugate(const TTD_CARRAY(TWO_TIMES_RANK)& array)
{
  BOOST_STATIC_ASSERT( tmptools::IsEvenAssert<TWO_TIMES_RANK>::value );
  TTD_CARRAY(TWO_TIMES_RANK) res(conj(array));
  details::helpHermitianConjugate(res);
  return res;
}


namespace dodirect {

using namespace linalg;

template<bool> void doDirect(CMatrix&, const CVector&, const CVector&);


struct Mul : boost::mpl::true_  {};
struct Add : boost::mpl::false_ {};


} // dodirect



template<int RANK1, int RANK2, bool MULT>
const TTD_CARRAY(RANK1+RANK2)
doDirect(const TTD_CARRAY(RANK1)& array1, const TTD_CARRAY(RANK2)& array2, boost::mpl::bool_<MULT>)
{
  using namespace linalg;
  TTD_CARRAY(RANK1+RANK2) res(concatenateTinies(array1.shape(),array2.shape()));
  if (res.data()) {
    CVector r1array1(unaryArray(array1));
    CVector r1array2(unaryArray(array2));
    CMatrix r2res(res.data(),shape(array1.size(),array2.size()),neverDeleteData);
    
    dodirect::doDirect<MULT>(r2res,r1array1,r1array2);
  }
  return res;

}



} // blitzplusplus

#endif // UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED
