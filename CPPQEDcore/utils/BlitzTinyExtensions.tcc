// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_TCC_INCLUDED

#include "BlitzTinyExtensions.h"

#include "Conversions.h"
#include "TMP_Tools.h"

#include <boost/mpl/for_each.hpp>



namespace blitzplusplus {


using blitz::TinyVector;


namespace details {


template<typename V1, typename V2, int OFFSET=0>
class TinyAssigner
{
public:
  TinyAssigner(V1& v1, const V2& v2) : v1_(v1), v2_(v2) {}

  template<typename T>
  void operator()(T)
  {
    typedef typename V1::T_numtype V1_numtype;
    typedef typename V2::T_numtype V2_numtype;

    v1_(T::value+OFFSET)=cpputils::Converter<V1_numtype,V2_numtype>::convert(v2_(T::value));
  }

private:
  V1& v1_;
  const V2& v2_;
  
};


} // details



template<typename T1, typename T2, int RANK1, int RANK2>
TinyVector<T1,RANK1+RANK2>
concatenateTinies(const TinyVector<T1,RANK1>& op1, const TinyVector<T2,RANK2>& op2)
{
  using boost::mpl::for_each;

  typedef TinyVector<T1,RANK1> OP1;
  typedef TinyVector<T2,RANK2> OP2;
  typedef TinyVector<T1,RANK1+RANK2> RES;

  RES res;

  details::TinyAssigner<RES,OP1      > assigner1(res,op1);
  details::TinyAssigner<RES,OP2,RANK1> assigner2(res,op2);

  for_each<tmptools::Ordinals<RANK1> >(assigner1);
  for_each<tmptools::Ordinals<RANK2> >(assigner2);

  return res;
    
}


namespace details {


#ifndef   NDEBUG

template<typename T, int RANK>
class TinyChecker
{
public:
  TinyChecker(const TinyVector<T,2*RANK>& v) : v_(v) {}

  template<typename U>
  void operator()(U)
  {
    if (v_(U::value)!=v_(U::value+RANK)) throw HalfCutTinyException();
  }

private:
  const TinyVector<T,2*RANK>& v_;
  
};

#endif // NDEBUG


} // details


template<typename T, int TWO_TIMES_RANK> 
TinyVector<T,TWO_TIMES_RANK/2> 
halfCutTiny(const TinyVector<T,TWO_TIMES_RANK>& tiny)
{
  using boost::mpl::for_each;

  static const int RANK=tmptools::AssertEvenAndDivideBy2<TWO_TIMES_RANK>::value;

#ifndef   NDEBUG
  details::TinyChecker<T,RANK> checker(tiny);
  boost::mpl::for_each<tmptools::Ordinals<RANK> >(checker);
#endif // NDEBUG

  typedef TinyVector<T,TWO_TIMES_RANK> OP ;
  typedef TinyVector<T,          RANK> RES;

  RES res;

  details::TinyAssigner<RES,OP> assigner(res,tiny);

  for_each<tmptools::Ordinals<RANK> >(assigner);

  return res;

}


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_TCC_INCLUDED
