// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::ElementLiouvillean}
#ifndef STRUCTURE_ELEMENTLIOUVILLEAN_TCC_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEAN_TCC_INCLUDED

#include "ElementLiouvillean.h"

#include <boost/range/algorithm/transform.hpp>

#include <boost/bind.hpp>


namespace structure { namespace details {

template<int RANK, int NLINDBLADS>
inline
const boost::function<double(typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,false>::JumpRateStrategy)>
convert( NoTime  , const typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,false>::LazyDensityOperator& matrix)
{
  return boost::bind(&ElementLiouvilleanStrategies<RANK,NLINDBLADS,false>::JumpRateStrategy::operator(),_1  ,boost::cref(matrix));
}


template<int RANK, int NLINDBLADS>
inline
const boost::function<double(typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,true >::JumpRateStrategy)>
convert(OneTime t, const typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,true >::LazyDensityOperator& matrix)
{
  return boost::bind(&ElementLiouvilleanStrategies<RANK,NLINDBLADS,true >::JumpRateStrategy::operator(),_1,t,boost::cref(matrix));
}


template<int RANK, int NLINDBLADS>
inline void performJump( NoTime  , typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,false>::StateVectorLow& psi, typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,false>::JumpStrategy jump)
{
  jump(  psi);
}
  

template<int RANK, int NLINDBLADS>
inline void performJump(OneTime t, typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,true >::StateVectorLow& psi, typename ElementLiouvilleanStrategies<RANK,NLINDBLADS,true >::JumpStrategy jump)
{
  jump(t,psi);
}


} } // structure::details


template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
auto structure::ElementLiouvilleanStrategies<RANK,NLINDBLADS,IS_TIME_DEPENDENT>::average_v(Time t, const LazyDensityOperator& matrix) const -> const Rates
{
  Rates rates(NLINDBLADS); // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

  boost::transform(jumpRates_,rates.begin(),details::convert<RANK,NLINDBLADS>(t,matrix));
  
  return rates;
  
}


template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
void structure::ElementLiouvilleanStrategies<RANK,NLINDBLADS,IS_TIME_DEPENDENT>::actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const
{
  details::performJump<RANK,NLINDBLADS>(t,psi,jumps_(lindbladNo));
}


#endif // STRUCTURE_ELEMENTLIOUVILLEAN_TCC_INCLUDED


