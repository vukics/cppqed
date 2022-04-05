// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "TMP_Tools.h"

#include <array>
#include <concepts>
#include <span>
#include <vector>


namespace cppqedutils {


template <typename T, unsigned RANK>
class MultiArrayView
{
public:
  /*
  template <typename... Indices>
  requires ( sizeof...(Indices)==RANK && ( ... && std::is_convertible_v<Indices,size_t> ) )
  T& operator() (Indices... i);
  */
  
  T& operator() (std::convertible_to<size_t> auto ... i) requires (sizeof...(i)==RANK)
  {
#ifndef   NDEBUG
    hana::for_each(tmptools::ordinals<RANK>,[](const auto&) {});
#endif // NDEBUG
  }
  
private:
  size_t offset_;
  
  std::array<size_t,RANK> extents_, strides_;
  
  std::span<T> dataView_;
  
};


template <typename T, unsigned RANK>
class MultiArray : public MultiArrayView<T,RANK>
{
private:
  std::vector<T> data_;
  
};


} // cppqedutils
