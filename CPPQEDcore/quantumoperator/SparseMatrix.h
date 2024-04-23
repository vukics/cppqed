// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"

#include <list>


namespace quantumoperator {

using namespace ::quantumdata;


struct SparseMatrix
{
  using Elements = std::list<std::tuple<size_t,size_t,dcomp>>;
  // TODO: std::unordered_map<std::pair<size_t,size_t>,dcomp>; would be ideal, but the pair doesn’t seem to be hashable

  Elements elements;

  /// Applying as a Hamiltonian
  template <size_t RANK>
  void operator () (StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt) const
  {
    this->operator()(psi.dataView,dpsidt.dataView);
  }

  void operator () (std::span<const dcomp> psi, std::span<dcomp> dpsidt) const;

};


} // quantumoperator
