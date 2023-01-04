// Copyright András Vukics 2022–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "MultiArray.h"


namespace quantumoperator {
  

/// `ListOfOffsets` should be a `hana::tuple` of `hana::tuple`s of `size_t`s
template <typename ListOfOffsets>
class MultiDiagonal
{
public:
  static constexpr size_t N_RANK=hana::size(ListOfOffsets{});
  
};
  
  
} // quantumoperator
