// -*- C++ -*-

#define n BOOST_PP_ITERATION()

#define TOA_print1(z,m,unused) (*this)(m).reference(x##m);
#define TOA_print2(z,m,unused) (*this)(m).reference(x##m.copy());


TinyOfArrays(TOA_ShallowCopy, BOOST_PP_ENUM_PARAMS(n,const T_numtype& x) )
{
  BOOST_PP_REPEAT(n,TOA_print1,~);
}

TinyOfArrays(TOA_DeepCopy   , BOOST_PP_ENUM_PARAMS(n,const T_numtype& x) )
{
  BOOST_PP_REPEAT(n,TOA_print2,~);
}


#undef TOA_print2
#undef TOA_print1

#undef n
