// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"

#include "MultiDiagonal.h"


namespace mode {


using MultiDiagonal = ::quantumoperator::MultiDiagonal<1> ;


MultiDiagonal aOp(size_t cutoff);

MultiDiagonal aDagOp(size_t cutoff);

MultiDiagonal nOp(size_t cutoff);

} // mode
