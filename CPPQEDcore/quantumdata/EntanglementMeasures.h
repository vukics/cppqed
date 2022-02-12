// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Declares the function calculating the negativity of a partially transposed density operator}
#ifndef   CPPQEDCORE_QUANTUMDATA_NEGPT_H_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_NEGPT_H_INCLUDED

#include "DensityOperator.h"

#ifdef EIGEN3_FOUND

#include "SliceIterator.h"
#include "TMP_Tools.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 


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
template<int RANK, typename V>
double negPT(const DensityOperator<RANK>& rho, V)
{
  using namespace boost::mpl;
  namespace mpl=boost::mpl;
  using namespace tmptools;

  using ExtendV = typename
    fold<Range<RANK,RANK>,
        typename fold<Ordinals<RANK>,
                      vector_c<int>,
                      push_back<mpl::_1,
                                if_<numerical_contains<V,mpl::_2>,
                                    plus<mpl::_2,int_<RANK> >,
                                    mpl::_2
                                    >
                                >
                      >::type,
        push_back<mpl::_1,
                  if_<numerical_contains<V,minus<mpl::_2,int_<RANK> > >,
                      minus<mpl::_2,int_<RANK> >,
                      mpl::_2
                      >
                  >
        >::type;

  DensityOperatorLow<RANK> rhoShallowPT(rho.getArray());

  cppqedutils::sliceiterator::Transposer<CArray,2*RANK,ExtendV>::_(rhoShallowPT);

  DensityOperatorLow<RANK> rhoDeepPT(rhoShallowPT.shape()); rhoDeepPT=rhoShallowPT;
  
  auto ev=Eigen::ComplexEigenSolver<Eigen::MatrixX<dcomp>>{
    Eigen::Map<Eigen::MatrixX<dcomp>>{rhoDeepPT.data(),long(rho.getTotalDimension()),long(rho.getTotalDimension())},
    false}.eigenvalues();
  
  return std::accumulate(ev.begin(),ev.end(),0.,[&] (double v, dcomp e) {return v - (real(e)<0 ? real(e) : 0.) ;}) ;
    
}

template<int RANK>
inline
double negPT(const DensityOperator<RANK>&, tmptools::V_Empty)
{
  return 0;
}


} // quantumdata 


#else  // EIGEN3_FOUND

namespace quantumdata {

/** If FLENS is not used, a dummy definition of negPT is provided. */
template<int RANK, typename V>
inline
std::string negPT(const DensityOperator<RANK>&, V)
{
  return "n/a";
}


} // quantumdata 


#endif // EIGEN3_FOUND


#endif // CPPQEDCORE_QUANTUMDATA_NEGPT_H_INCLUDED
