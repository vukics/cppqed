// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Range algorithms not yet in STL in C++20}
#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>


namespace cppqedutils {


/// std::array`s of different size cannot be std::views::join`ed because the size is part of the type
/** this solution based on C++17 fold expressions comes from [here](http://stackoverflow.com/a/42774523/1171157) */
template <typename Type, std::size_t... sizes>
constexpr auto concatenate(const std::array<Type, sizes>&... arrays)
{
  std::array<Type, (sizes + ...)> result;
  std::size_t index{};

  ((std::copy_n(arrays.begin(), sizes, result.begin() + index), index += sizes), ...);

  return result;
}


} // cppqedutils
