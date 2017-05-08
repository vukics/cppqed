// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Tridiagonal.h"

#define RANK 7


typedef quantumoperator::Tridiagonal<RANK>::StateVectorLow StateVectorLow;

template void quantumoperator::Tridiagonal<RANK>::apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const;

