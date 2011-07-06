// -*- C++ -*-
#ifndef   TRIDIAGONAL_HAMILTONIAN_IMPL_INCLUDED
#define   TRIDIAGONAL_HAMILTONIAN_IMPL_INCLUDED

#include "Algorithm.h"

#include "Range.h"

#include <boost/bind.hpp>


namespace structure {


template<int RANK>
void 
details::TDH_True<RANK>::addContribution(double tMinusIntPic0, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,bind(&Tridiagonal::propagate,_1,tMinusIntPic0)));
}


template<int RANK>
void 
details::TDH_False<RANK>::addContribution(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,_1));
}



} // structure

#endif // TRIDIAGONAL_HAMILTONIAN_IMPL_INCLUDED
