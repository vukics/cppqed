#include "Tridiagonal.h"

#define RANK 5


typedef typename quantumoperator::Tridiagonal<RANK>::StateVectorLow StateVectorLow;

template void quantumoperator::Tridiagonal<RANK>::apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const;

