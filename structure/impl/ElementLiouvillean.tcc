// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::ElementLiouvillean}
#ifndef STRUCTURE_IMPL_ELEMENTLIOUVILLEAN_TCC_INCLUDED
#define STRUCTURE_IMPL_ELEMENTLIOUVILLEAN_TCC_INCLUDED

#include "ElementLiouvillean.h"

#include "Range.h"


namespace structure { namespace details {

template<int RANK, int NOJ>
inline
const boost::function<double(typename ElementLiouvilleanStrategies<RANK,NOJ,false>::JumpRateStrategy)>
convert( NoTime  , const typename ElementLiouvilleanStrategies<RANK,NOJ,false>::LazyDensityOperator& matrix)
{
  return boost::bind(&ElementLiouvilleanStrategies<RANK,NOJ,false>::JumpRateStrategy::operator(),_1  ,boost::cref(matrix));
}


template<int RANK, int NOJ>
inline
const boost::function<double(typename ElementLiouvilleanStrategies<RANK,NOJ,true >::JumpRateStrategy)>
convert(OneTime t, const typename ElementLiouvilleanStrategies<RANK,NOJ,true >::LazyDensityOperator& matrix)
{
  return boost::bind(&ElementLiouvilleanStrategies<RANK,NOJ,true >::JumpRateStrategy::operator(),_1,t,boost::cref(matrix));
}


template<int RANK, int NOJ>
inline void performJump( NoTime  , typename ElementLiouvilleanStrategies<RANK,NOJ,false>::StateVectorLow& psi, typename ElementLiouvilleanStrategies<RANK,NOJ,false>::JumpStrategy jump)
{
  jump(  psi);
}
  

template<int RANK, int NOJ>
inline void performJump(OneTime t, typename ElementLiouvilleanStrategies<RANK,NOJ,true >::StateVectorLow& psi, typename ElementLiouvilleanStrategies<RANK,NOJ,true >::JumpStrategy jump)
{
  jump(t,psi);
}


} } // structure::details


template<int RANK, int NOJ, bool IS_TIME_DEPENDENT>
auto structure::ElementLiouvilleanStrategies<RANK,NOJ,IS_TIME_DEPENDENT>::average_v(Time t, const LazyDensityOperator& matrix) const -> const Rates
{
  Rates rates(NOJ); // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

  boost::transform(jumpRates_,rates.begin(),details::convert<RANK,NOJ>(t,matrix));
  
  return rates;
  
}


template<int RANK, int NOJ, bool IS_TIME_DEPENDENT>
void structure::ElementLiouvilleanStrategies<RANK,NOJ,IS_TIME_DEPENDENT>::actWithJ_v(Time t, StateVectorLow& psi, size_t jumpNo) const
{
  details::performJump<RANK,NOJ>(t,psi,jumps_(jumpNo));
}


#endif // STRUCTURE_IMPL_ELEMENTLIOUVILLEAN_TCC_INCLUDED


