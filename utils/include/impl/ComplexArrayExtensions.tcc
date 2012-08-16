// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_COMPLEXARRAYEXTENSIONS_TCC_INCLUDED

#include "ComplexArrayExtensions.h"

#include "impl/BlitzArrayExtensions.tcc"
#include "impl/BlitzTinyExtensions.tcc"


namespace blitzplusplus {


namespace details {

template<int RANK>
void helpHermitianConjugate(TTD_CARRAY(2*RANK)&);

} // details

template<int TWO_TIMES_RANK>
inline
void
hermitianConjugateSelf(TTD_CARRAY(TWO_TIMES_RANK)& array) 
{
  static const int RANK=tmptools::IsEvenAssert<TWO_TIMES_RANK>::value;
  array=conj(array); 
  details::helpHermitianConjugate<RANK>(array);
}


template<int TWO_TIMES_RANK>
inline
const TTD_CARRAY(TWO_TIMES_RANK)
hermitianConjugate(const TTD_CARRAY(TWO_TIMES_RANK)& array)
// NEED_TO_UNDERSTAND Funnily enough, if it is declared in such a way:
// hermitianConjugate(const typename CAMF<TWO_TIMES_RANK>::type&)
// then the compiler cannot deduce the argument any more
{
  static const int RANK=tmptools::IsEvenAssert<TWO_TIMES_RANK>::value;
  TTD_CARRAY(TWO_TIMES_RANK) res(conj(array));
  details::helpHermitianConjugate<RANK>(res);
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
