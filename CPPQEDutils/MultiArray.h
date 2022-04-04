// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include <array>
#include <vector>


using cvector=std::vector<const double>;


template <typename T, unsigned RANK, template <typename> class STORAGE>
class MultiArray
{
public:
  template <typename... Indices>
  requires ( sizeof...(Indices)==RANK ) // + all indices should be (convertible to) size_t
  T& operator() (Indices... i);
  
private:
  size_t offset_;
  
  std::array<size_t,RANK> extents_, strides_;
  
};
