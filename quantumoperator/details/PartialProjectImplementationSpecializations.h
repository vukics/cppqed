// -*- C++ -*-

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
