// -*- C++ -*-
#ifndef _NEGPT_H
#define _NEGPT_H

#include "DensityOperatorFwd.h"

#include "TMP_Tools.h"


namespace quantumdata {

// Calculates the negativity of the partial transpose of the density
// operator of an arbitrarily complex system. Of course it should be
// regarded as a bipartite system, so that a subsystem has to be
// specified to be one party of the bipartite.

// V should be a compile-time vector, tipically a tmptools::vector
// specifying the subsystem

template<int RANK, typename V>
double negPT(const DensityOperator<RANK>&, V=V());

template<int RANK>
inline
double negPT(const DensityOperator<RANK>&, tmptools::V0)
{
  return 0;
}


} // quantumdata 

#include "impl/NegPT.tcc"

#endif // _NEGPT_H
