// -*- C++ -*-
#ifndef   STRUCTURE_IMPL_TRIDIAGONALHAMILTONIAN_TCC_INCLUDED
#define   STRUCTURE_IMPL_TRIDIAGONALHAMILTONIAN_TCC_INCLUDED

#include "TridiagonalHamiltonian.h"

#include "Algorithm.h"
#include "Range.h"
#include "impl/Tridiagonal.tcc"

#include <boost/bind.hpp>


namespace structure {


template<int RANK>
void 
details::TDH_True<RANK>::addContribution_v(double tMinusIntPic0, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,bind(&Tridiagonal::propagate,_1,tMinusIntPic0)));
}


template<int RANK>
void 
details::TDH_False<RANK>::addContribution_v(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,_1));
}



} // structure

#endif // STRUCTURE_IMPL_TRIDIAGONALHAMILTONIAN_TCC_INCLUDED
