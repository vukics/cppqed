// -*- C++ -*-

#define N BOOST_PP_ITERATION()

template<BOOST_PP_ENUM_PARAMS(N,typename T)>
typename SliceInfo<T_numtype,BOOST_PP_ENUM_PARAMS(N,T)>::T_slice
operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,T,r)) const
{
  typedef typename SliceInfo<T_numtype,BOOST_PP_ENUM_PARAMS(N,T)>::T_slice slice;
  return slice(noConst(), BOOST_PP_ENUM_PARAMS(N,r) BOOST_PP_ENUM_TRAILING(BOOST_PP_SUB(BLITZ_ARRAY_LARGEST_RANK,N),DEFAULT_print,~) );
}
