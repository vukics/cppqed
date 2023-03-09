// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include <type_traits>

namespace cppqedutils {

/// The overload pattern for std::visit of std::variants, cf. https://www.cppstories.com/2018/09/visit-variants/
template<class... Ts> struct overload : Ts... { using Ts::operator()...; };


/// An empty base class with a constructor of arbitrary signature
struct Empty { template <typename... T> Empty(T&&... ) {} };


template <typename T>
constexpr bool passByValue_v=false;


template <typename State> struct ReferenceMF : std::type_identity<std::add_lvalue_reference_t<State>> {};

template <typename State> using Reference = typename ReferenceMF<State>::type;

template <typename State> struct ConstReferenceMF : std::type_identity<std::add_lvalue_reference_t<std::add_const_t<State>>> {};

template <typename State> using ConstReference = typename ConstReferenceMF<State>::type;

} // cppqedutils
