// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the hierarchical partial specializations of structure::ElementLiouvillean}
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_TCC_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_TCC_INCLUDED

#include "ElementLiouvillean.h"

#include <boost/range/algorithm/transform.hpp>

#include <boost/bind.hpp>


namespace structure { namespace details {

template<int RANK>
inline
const boost::function<double(typename ElementLiouvilleanStrategiesBase<RANK,false>::JumpRateStrategy)>
convert( NoTime  , const typename ElementLiouvilleanStrategiesBase<RANK,false>::LazyDensityOperator& matrix)
{
  return boost::bind(&ElementLiouvilleanStrategiesBase<RANK,false>::JumpRateStrategy::operator(),_1  ,boost::cref(matrix));
}


template<int RANK>
inline
const boost::function<double(typename ElementLiouvilleanStrategiesBase<RANK,true >::JumpRateStrategy)>
convert(OneTime t, const typename ElementLiouvilleanStrategiesBase<RANK,true >::LazyDensityOperator& matrix)
{
  return boost::bind(&ElementLiouvilleanStrategiesBase<RANK,true >::JumpRateStrategy::operator(),_1,t,boost::cref(matrix));
}


template<int RANK>
inline void performJump( NoTime  , typename ElementLiouvilleanStrategiesBase<RANK,false>::StateVectorLow& psi, typename ElementLiouvilleanStrategiesBase<RANK,false>::JumpStrategy jump)
{
  jump(  psi);
}
  

template<int RANK>
inline void performJump(OneTime t, typename ElementLiouvilleanStrategiesBase<RANK,true >::StateVectorLow& psi, typename ElementLiouvilleanStrategiesBase<RANK,true >::JumpStrategy jump)
{
  jump(t,psi);
}


template<int RANK>
inline void performSuperoperation(NoTime,
                                  const typename ElementLiouvilleanStrategiesBase<RANK,false>::DensityOperatorLow& rho,
                                  typename ElementLiouvilleanStrategiesBase<RANK,false>::DensityOperatorLow& drhodt,
                                  typename ElementLiouvilleanStrategiesBase<RANK,false>::SuperoperatorStrategy f
                                 )
{f(rho,drhodt);}


template<int RANK>
inline void performSuperoperation(OneTime t,
                                  const typename ElementLiouvilleanStrategiesBase<RANK,true>::DensityOperatorLow& rho,
                                  typename ElementLiouvilleanStrategiesBase<RANK,true>::DensityOperatorLow& drhodt,
                                  typename ElementLiouvilleanStrategiesBase<RANK,true>::SuperoperatorStrategy f
                                 )
{f(t,rho,drhodt);}


} } // structure::details


template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
auto structure::ElementLiouvilleanStrategies<RANK,NLINDBLADS,IS_TIME_DEPENDENT>::rates_v(Time t, const LazyDensityOperator& matrix) const -> const Rates
{
  Rates rates(NLINDBLADS); // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

  boost::transform(jumpRates_,rates.begin(),details::convert<RANK>(t,matrix));
  
  return rates;
  
}


template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
void structure::ElementLiouvilleanStrategies<RANK,NLINDBLADS,IS_TIME_DEPENDENT>::actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const
{
  details::performJump<RANK>(t,psi,jumps_(lindbladNo));
}


template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
void structure::ElementLiouvilleanStrategies<RANK,NLINDBLADS,IS_TIME_DEPENDENT>::actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const
{
  if (superoperatorStrategies_(lindbladNo).empty()) throw SuperoperatorNotImplementedException(lindbladNo);
  details::performSuperoperation<RANK>(t,rho,drhodt,superoperatorStrategies_(lindbladNo));
}

#endif // CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_TCC_INCLUDED


