// -*- C++ -*-
#ifndef   NEGPT_H_INCLUDED
#define   NEGPT_H_INCLUDED

#include "DensityOperator.h"

#include "BlitzArraySliceIterator.h"
#include "TMP_Tools.h"



#ifndef   DO_NOT_USE_FLENS


namespace quantumdata {

// Calculates the negativity of the partial transpose of the density
// operator of an arbitrarily complex system. Of course it should be
// regarded as a bipartite system, so that a subsystem has to be
// specified to be one party of the bipartite.

// V should be a compile-time vector, tipically a tmptools::vector
// specifying the subsystem

template<int RANK, typename V>
double negPT(const DensityOperator<RANK>&, V);

template<int RANK>
inline
double negPT(const DensityOperator<RANK>&, tmptools::V_Empty)
{
  return 0;
}


} // quantumdata 


#include "impl/NegPT.tcc"


#else  // DO_NOT_USE_FLENS

namespace quantumdata {

// If FLENS is not used, a dummy definition of negPT is provided

template<int RANK, typename V>
double negPT(const DensityOperator<RANK>&, V)
{
  return 0;
}


} // quantumdata 


#endif // DO_NOT_USE_FLENS


#endif // NEGPT_H_INCLUDED
