// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ComplexExtensions.h"

#include <type_traits>
#include <valarray>

namespace cppqedutils {

/// The overload pattern for std::visit of std::variants, cf. https://www.cppstories.com/2018/09/visit-variants/
template<class... Ts> struct overload : Ts... { using Ts::operator()...; };


/// An empty base class with a constructor of arbitrary signature
struct Empty { template <typename... T> Empty(T&&... ) {} };


template <typename T>
constexpr bool passByValue_v=false;


/// metafunctions are needed since template aliases cannot be partially specialized
template <typename State> struct ReferenceMF : std::type_identity<std::add_lvalue_reference_t<State>> {};

template <typename State> using Reference = typename ReferenceMF<State>::type;

template <typename State> struct ConstReferenceMF : std::type_identity<std::add_lvalue_reference_t<std::add_const_t<State>>> {};

template <typename State> using ConstReference = typename ConstReferenceMF<State>::type;


/// either a valarray of dcomp or valarray of such (recursive definition)
/** Workaround for defining concept recursively from [this Q&A](https://stackoverflow.com/questions/56741456/how-to-define-a-recursive-concept) */
namespace traits {

template <typename T>
constexpr bool tdp = false;

template <>
constexpr bool tdp<std::valarray<dcomp>> = true;

template <typename T>
constexpr bool tdp<std::valarray<T>> = tdp<T>;

} // traits

template <typename T> concept temporal_data_point = traits::tdp<T>;


static_assert(temporal_data_point<std::valarray<dcomp>>);
static_assert(temporal_data_point<std::valarray<std::valarray<dcomp>>>);
static_assert(temporal_data_point<std::valarray<std::valarray<std::valarray<dcomp>>>>);
static_assert(!temporal_data_point<std::valarray<std::valarray<std::valarray<int>>>>);


} // cppqedutils



#include <boost/hana.hpp>

namespace hana=boost::hana;

template <typename S> concept hana_sequence = hana::Sequence<S>::value;


#include <boost/json.hpp>

namespace cppqedutils {

using LogTree = boost::json::object;

template <typename T> concept intro_logger = requires ( const T& t ) {
  { logIntro(t) } -> std::convertible_to<LogTree>; };

template <typename T> concept outro_logger = requires ( const T& t ) {
  { logOutro(t) } -> std::convertible_to<LogTree>; };

template <typename T> concept logger = intro_logger<T> && outro_logger<T>;

} // cppqedutils
