#include "Tridiagonal.h"

#define RANK 7


typedef quantumoperator::Tridiagonal<RANK>::StateVectorLow StateVectorLow;

template void quantumoperator::Tridiagonal<RANK>::apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const;

