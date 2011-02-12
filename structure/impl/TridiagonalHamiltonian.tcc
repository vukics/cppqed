// -*- C++ -*-
#ifndef   TRIDIAGONAL_HAMILTONIAN_IMPL_INCLUDED
#define   TRIDIAGONAL_HAMILTONIAN_IMPL_INCLUDED

#include "Algorithm.h"

#include "Range.h"

#include <boost/iterator/transform_iterator.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace structure {


template<int RANK>
void 
details::TDH_True<RANK>::addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  cpputils::for_each(hOverIs_,freqss_.begin(),bind(quantumoperator::apply<RANK>,psi,dpsidt,bind(&Tridiagonal::propagate,_1,_2,t-tIntPic0)));
}


template<int RANK>
void 
details::TDH_False<RANK>::addContribution(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  boost::for_each(hOverIs_,bind(quantumoperator::apply<RANK>,psi,dpsidt,_1));
}


namespace details {


template<int RANK, int GETN, typename TS>
const TS
extractor(const typename TridiagonalHamiltonian<RANK,true>::TridiagonalIPs& tridiagonalIPs)
{
  using namespace std  ;
  using namespace boost;

  typedef typename TridiagonalHamiltonian<RANK,true>::TridiagonalIP TridiagonalIP;
  
  typedef typename TS::value_type VT;

  TS res;

  function<const VT&(const TridiagonalIP&)> transformation(mem_fn<const VT&>(&TridiagonalIP::template get<GETN>)
							   /*bind(&TridiagonalIP::template get<GETN>,_1)*/);

  copy(make_transform_iterator(tridiagonalIPs.begin(),transformation),make_transform_iterator(tridiagonalIPs.end(),transformation),
       back_inserter(res));

  return res;

}



} // details


template<int RANK, bool IS_TD>
TridiagonalHamiltonian<RANK,IS_TD>::TridiagonalHamiltonian(const typename TridiagonalHamiltonian<RANK,IS_TD>::TridiagonalIPs& tridiagonalIPs)
  : Base(details::extractor<RANK,0,typename TridiagonalHamiltonian<RANK,true>::Tridiagonals>(tridiagonalIPs),
	 details::extractor<RANK,1,typename TridiagonalHamiltonian<RANK,true>::Frequenciess>(tridiagonalIPs))
{
}


} // structure

#endif // TRIDIAGONAL_HAMILTONIAN_IMPL_INCLUDED
