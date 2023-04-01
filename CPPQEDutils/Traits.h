// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ComplexExtensions.h"

#include <boost/hana.hpp>

#include <type_traits>


namespace hana=boost::hana;

template <typename S> concept hana_sequence = hana::Sequence<S>::value;


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


/// hana::tuple of double/dcomp or such hana::tuples (recursive definition)
/** Workaround for defining concept recursively from [this Q&A](https://stackoverflow.com/questions/56741456/how-to-define-a-recursive-concept) */
namespace traits {

template <typename T> constexpr bool tdp = false;

template <> constexpr bool tdp<double> = true;

template <> constexpr bool tdp<dcomp> = true;

template <hana_sequence S> constexpr bool tdp<S> = !!hana::all_of(
  decltype(hana::transform(std::declval<S>(), hana::typeid_)){},
  []<class T>(T) { return tdp<typename T::type>; });

} // traits

template <typename T> concept temporal_data_point = traits::tdp<T>;


static_assert(temporal_data_point<double>);
static_assert(temporal_data_point<dcomp>);
static_assert(temporal_data_point<decltype(hana::make_tuple(1.,dcomp{2.,-1.}))>);
static_assert(temporal_data_point<decltype(hana::make_tuple(1.,dcomp{2.,-1.}),hana::make_tuple(1.,dcomp{2.,-1.}))>);
static_assert(!temporal_data_point<int>);


} // cppqedutils


#include <json.hpp>

namespace cppqedutils {

using LogTree = nlohmann::json ;

template <typename T> concept intro_logger = requires ( const T& t ) {
  { logIntro(t) } -> std::convertible_to<LogTree>; };

template <typename T> concept outro_logger = requires ( const T& t ) {
  { logOutro(t) } -> std::convertible_to<LogTree>; };

template <typename T> concept logger = intro_logger<T> && outro_logger<T>;

template <typename H>
concept labelled = requires (H&& h) { { label(h) } -> std::convertible_to<LogTree>; } || requires (H&& h) { { h.label } -> std::convertible_to<LogTree>; } ;


} // cppqedutils



/// cf. [this discussion](https://github.com/boostorg/hana/issues/317)
/** TODO: Unfortunately, this doesn’t seem to be enough for use in std::apply fold expressions */
// namespace std
// {
//     template<std::size_t n, typename... Types>
//     struct tuple_element<n, boost::hana::tuple<Types...>>
//     {
//         using type = typename decltype(+boost::hana::tuple_t<Types...>[boost::hana::size_c<n>])::type;
//     };
//
//     template<typename... Types>
//     struct tuple_size<boost::hana::tuple<Types...>>:
//         public integral_constant<std::size_t, sizeof...(Types)>
//     {};
// }
//
//
// namespace boost
// {
//     namespace hana
//     {
//         template<std::size_t n, typename... Types>
//         constexpr decltype(auto) get(hana::tuple<Types...>& t)
//         {
//             return t[hana::size_c<n>];
//         }
//
//         template<std::size_t n, typename... Types>
//         constexpr decltype(auto) get(const hana::tuple<Types...>& t)
//         {
//             return t[hana::size_c<n>];
//         }
//
//         template<std::size_t n, typename... Types>
//         constexpr decltype(auto) get(hana::tuple<Types...>&& t)
//         {
//             return static_cast<hana::tuple<Types...>&&>(t)[hana::size_c<n>];
//         }
//
//         template<std::size_t n, typename... Types>
//         constexpr decltype(auto) get(const hana::tuple<Types...>&& t)
//         {
//             return static_cast<const hana::tuple<Types...>&&>(t)[hana::size_c<n>];
//         }
//     }
// }
