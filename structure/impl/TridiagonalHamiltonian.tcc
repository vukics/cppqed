// -*- C++ -*-
#ifndef   STRUCTURE_IMPL_TRIDIAGONALHAMILTONIAN_TCC_INCLUDED
#define   STRUCTURE_IMPL_TRIDIAGONALHAMILTONIAN_TCC_INCLUDED

#include "TridiagonalHamiltonian.h"

#include "Algorithm.h"
#include "Range.h"
#include "impl/Tridiagonal.tcc"

#include <boost/utility/enable_if.hpp>
#include <boost/bind.hpp>


namespace structure {


template<int RANK, bool IS_TIME_DEPENDENT>
void 
details::TDH_Base<RANK,IS_TIME_DEPENDENT>::addContribution_v(OneTime t, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  typename boost::enable_if_c< IS_TIME_DEPENDENT>::type(); ///< Ensuring that the function can be compiled only in the time-dependent case
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,bind(&Tridiagonal::propagate,_1,double(t))));
}


template<int RANK, bool IS_TIME_DEPENDENT>
void 
details::TDH_Base<RANK,IS_TIME_DEPENDENT>::addContribution_v(NoTime, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  typename boost::enable_if_c<!IS_TIME_DEPENDENT>::type(); ///< Ensuring that the function can be compiled only in the time-independent case
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,_1));
}



} // structure

#endif // STRUCTURE_IMPL_TRIDIAGONALHAMILTONIAN_TCC_INCLUDED
