// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "DensityOperator.h"

#include <Eigen/Eigenvalues> 

#include <numeric>

namespace quantumdata {


  
/// Calculates the negativity of the partial transpose of the density operator of an arbitrarily complex system
/**
 * \see \cite vidal02
 *
 * The system should of course be regarded as a bipartite system, so that a subsystem has to be specified to be 
 * one part of the bipartite. The compile-time vector `V` specifies the subsystem.
 * 
 * The negativity is calculated as *the sum of the negative eigenvalues* of the partially transposed density operator.
 * 
 * \note This definition is equivalent to the original definition \cite vidal02 up to a sign,
 * because the partially transposed density operator's eigenvalues come in two sorts:
 * - solitary positive numbers (\f$a_i\f$) adding up to one, and
 * - pairs of opposite-sign numbers (\f$b_i\f$).
 * \note Hence
 * \f[\mathcal{N}(\rho)\equiv\frac{\norm{\rho^{\text{PT}}}_1-1}2=\frac{\sum_i\abs{\rho^{\text{PT}}_{ii}}-1}2=\frac{\sum_ia_i+2\sum_ib_i-1}2=\sum_ib_i.\f]\par
 * 
 * \tparamRANK
 * \tparam V a compile-time vector, tipically a tmptools::Vector, specifying the quantum numbers grouped into one part of the bipartite (cf. \ref specifyingsubsystems)
 * 
 */
template<auto subsystemAxes, size_t RANK>
double negPT(const DensityOperator<RANK>& rho)
{
  static constexpr auto extendedAxes = hana::fold (
    subsystemAxes ,
    cppqedutils::compileTimeOrdinals<2*RANK> ,
    [] (auto state, auto v) {state[v]=v+RANK; state[v+RANK]=v; return state;}
    );  

  DensityOperator<RANK> rhoDeepPT{rho.extents,
                                  [&] (auto& r) {r.mutableView.assignTo(cppqedutils::transpose<extendedAxes>(rho));}};
  
  auto ev=Eigen::ComplexEigenSolver<CMatrix>{Eigen::Map<CMatrix>{rhoDeepPT.data(),rho.getTotalDimension(),rho.getTotalDimension()},false}.eigenvalues();
  
  return (std::accumulate(ev.begin(),ev.end(),0.,[] (double v, dcomp e) {return v + std::abs(e) ;}) - 1.)/2. ;
  
}

template<size_t RANK> double negPT(const DensityOperator<RANK>&) {return 0;}


template<size_t RANK>
double entropy(const DensityOperator<RANK>& rho)
{
  auto ev=Eigen::SelfAdjointEigenSolver<CMatrix>{
    Eigen::Map<CMatrix>{const_cast<dcomp*>(rho.getArray().data()),long(rho.getTotalDimension()),long(rho.getTotalDimension())},Eigen::EigenvaluesOnly}.eigenvalues();
    
  return std::accumulate(ev.begin(),ev.end(),0.,[] (double v, double e) {return v - e * ( e>0. ? std::log(e) : 0. ) ;});
}


template<auto subsystemAxes, size_t RANK>
double mutualInformation(const DensityOperator<RANK>& rho)
{
  return entropy(reduce<V>(rho))
    +entropy(reduce<::tmptools::NegatedVector<RANK,V> >(rho))
    -entropy(rho);
}


template<int RANK, typename V>
double purityOfPartialTrace(const DensityOperator<RANK>& rho, V)
{
  return purity(reduce<V>(rho));
}



} // quantumdata 
