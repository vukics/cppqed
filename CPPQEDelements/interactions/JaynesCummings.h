// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Mode.h"


namespace jaynescummings {

using namespace ::quantumdata;


template <size_t i, size_t j>
auto hamiltonian(size_t cutoff, dcomp g)
{
  return [
    aOp = -g*::mode::aOp(cutoff),
    aDagOp = conj(g)*::mode::aDagOp(cutoff)
  ] (StateVectorConstView<2> psi, StateVectorView<2> dpsidt) {
    using ::cppqedutils::multiarray::_;
    aOp(psi(i,_),dpsidt(j,_));
    aDagOp(psi(j,_),dpsidt(i,_));
  };
}

} // jaynescummings



// JaynesCummings is a template which knows about the levels between which it steps @ compile time => it covers also the former MLCJ

// two version: one which treats the sigma operator “by hand”, and one with purely (binary) MultiDiagonal design
