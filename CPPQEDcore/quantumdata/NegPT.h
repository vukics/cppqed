// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Declares the function calculating the negativity of a partially transposed density operator}
#ifndef   CPPQEDCORE_QUANTUMDATA_NEGPT_H_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_NEGPT_H_INCLUDED

#include "DensityOperator.h"

#include "TMP_Tools.h"



#ifndef   DO_NOT_USE_FLENS


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
double negPT(const DensityOperator<RANK>&, V);

template<int RANK>
inline
double negPT(const DensityOperator<RANK>&, tmptools::V_Empty)
{
  return 0;
}


} // quantumdata 


#else  // DO_NOT_USE_FLENS

namespace quantumdata {

/** If FLENS is not used, a dummy definition of negPT is provided. */
template<int RANK, typename V>
inline
std::string negPT(const DensityOperator<RANK>&, V)
{
  return "n/a";
}


} // quantumdata 


#endif // DO_NOT_USE_FLENS


#endif // CPPQEDCORE_QUANTUMDATA_NEGPT_H_INCLUDED
