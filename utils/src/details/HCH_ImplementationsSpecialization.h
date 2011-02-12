// -*- C++ -*-

#define ITER   BOOST_PP_ITERATION()
#define ITERT2 BOOST_PP_MUL(2,ITER)

#define HCH_print(z,m,unused) m ,

namespace blitzplusplus {
namespace details {

template<> 
void 
helpHermitianConjugate<ITER>(TTD_CARRAY(ITERT2)& array)
{
  array.transposeSelf(BOOST_PP_REPEAT_FROM_TO(ITER,ITERT2,HCH_print,~) BOOST_PP_ENUM_PARAMS(ITER,));
}

}
}

#undef HCH_print

#undef ITERT2
#undef ITER
