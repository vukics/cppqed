// -*- C++ -*-
#ifndef   _BLITZ_TINY_EXTENSIONS_IMPL_H
#define   _BLITZ_TINY_EXTENSIONS_IMPL_H

#include "Conversions.h"
#include "TMP_Tools.h"

#include<boost/mpl/for_each.hpp>



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

    v1_(T::value+OFFSET)=cpputils::ConverterMF<V1_numtype,V2_numtype>::type::convert(v2_(T::value));
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

  for_each<typename tmptools::OrdinalMF<RANK1>::type>(assigner1);
  for_each<typename tmptools::OrdinalMF<RANK2>::type>(assigner2);

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
    if (v_(U::value)!=v_(U::value+RANK)) throw HalfCutTinyFishyException();
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

  static const int RANK=tmptools::IsEvenAssert<TWO_TIMES_RANK>::value;

#ifndef   NDEBUG
  details::TinyChecker<T,RANK> checker(tiny);
  boost::mpl::for_each<typename tmptools::OrdinalMF<RANK>::type>(checker);
#endif // NDEBUG

  typedef TinyVector<T,TWO_TIMES_RANK> OP ;
  typedef TinyVector<T,          RANK> RES;

  RES res;

  details::TinyAssigner<RES,OP> assigner(res,tiny);

  for_each<typename tmptools::OrdinalMF<RANK>::type>(assigner);

  return res;

}


} // blitzplusplus


#endif // _BLITZ_TINY_EXTENSIONS_IMPL_H
