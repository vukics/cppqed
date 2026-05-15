// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

/// @brief Perfect-forwarding macro. Expands to `std::forward<decltype(x)>(x)`, avoiding the need to repeat the type.
#define FWD(x) std::forward<decltype(x)>(x)

#include <boost/hana.hpp>
/// @brief Namespace alias for Boost.Hana, used for compile-time sequences and heterogeneous containers.
namespace hana=boost::hana;

#include <boost/json.hpp>

#include <boost/serialization/split_free.hpp>

#include <numeric>
#include <type_traits>


namespace cppqedutils {

/// @brief Concept satisfied by any Boost.Hana sequence (tuple-like heterogeneous container).
/// Used to constrain template parameters that must be traversable at compile time via `hana::for_each`, `hana::transform`, etc.
template <typename S> concept hana_sequence = hana::Sequence<S>::value;

/// @brief Namespace alias for Boost.JSON, used for structured logging and JSON serialization.
namespace json = boost::json ;
/// @brief Alias for a JSON object used as a structured log tree.
/// Passed through the trajectory and structure layers to accumulate human-readable, machine-parseable log data.
using LogTree = json::object ;

/// @brief Serializes any JSON-convertible value to a string via `json::value_from`.
/// @param v Any value for which `boost::json::value_from` is defined.
/// @return A JSON string representation of @p v.
std::string toStringJSON(auto&& v) {return json::serialize( json::value_from( FWD(v) ) ) ;}

/// @brief Helper to produce a `static_assert` that always fails, with a dependent type in the error message.
/// Use as `static_assert(always_false<T>::value, "...")` inside unreachable `if constexpr` branches.
template <typename T>
struct always_false : std::false_type {};

} // cppqedutils


namespace cppqedutils {

/// @brief The overload pattern for `std::visit` on `std::variant`: constructs a callable from a set of lambdas.
/// @see https://www.cppstories.com/2018/09/visit-variants/
template<class... Ts> struct overload : Ts... { using Ts::operator()...; };


/// @brief Variable template that is always `false`, for `static_assert` in `if constexpr` branches without a dependent type.
template<typename > inline constexpr bool always_false_v = false;


/// @brief An empty base class whose constructor accepts and ignores any arguments.
/// Useful as a no-op base in mixin hierarchies where a forwarding constructor call is required but has no effect.
struct Empty { Empty(auto&&... ) {} };


/// @brief Primary template for the pass-by-value trait.
/// Specializations are `true` for view types (e.g. `MultiArrayView`, `LazyDensityOperator`) and `false` for owning types (e.g. `MultiArray`).
/// The primary template returns `std::nullopt` to signal the absence of a specialization.
template <typename T>
constexpr auto passByValue_v=std::nullopt;


/// @brief Metafunction yielding the reference type of @p State. Defaults to `State&`.
/// Specialized for view types to return the view itself, since views are already reference-like.
/// @note Template aliases cannot be partially specialized, hence the metafunction class.
template <typename State> struct ReferenceMF : std::type_identity<std::add_lvalue_reference_t<State>> {};

/// @brief Convenience alias for `ReferenceMF<State>::type`.
template <typename State> using Reference = typename ReferenceMF<State>::type;

/// @brief Metafunction yielding the const-reference type of @p State. Defaults to `const State&`.
/// Specialized for view types to return a const view.
template <typename State> struct ConstReferenceMF : std::type_identity<std::add_lvalue_reference_t<std::add_const_t<State>>> {};

/// @brief Convenience alias for `ConstReferenceMF<State>::type`.
template <typename State> using ConstReference = typename ConstReferenceMF<State>::type;


/// @brief Concept for types that produce an intro log via ADL `logIntro`. Used in the trajectory layer.
template <typename T> concept intro_logger = requires ( const T& t ) {
  { logIntro(t) } -> std::convertible_to<LogTree>; };

/// @brief Concept for types that produce an outro log via ADL `logOutro` (e.g. final ODE step counts).
template <typename T> concept outro_logger = requires ( const T& t ) {
  { logOutro(t) } -> std::convertible_to<LogTree>; };

/// @brief Concept for types satisfying both `intro_logger` and `outro_logger`.
template <typename T> concept logger = intro_logger<T> && outro_logger<T>;


/// @brief Concept for types carrying a string label, either as a member `h.label` or via ADL `label(h)`.
/// Used in the structure and Liouvillian layers to name operators (e.g. Lindblad jump operators, expectation value columns).
/// @tparam OUT Type the label must be convertible to; defaults to `LogTree`.
template <typename H, typename OUT = LogTree>
concept labelled = requires (const H& h) { { h.label } -> std::convertible_to<OUT>; } || requires (const H& h) { { label(h) } -> std::convertible_to<OUT>; } ;


/// @brief Retrieves the label of a labelled object. Prefers the `h.label` member; falls back to ADL `label(h)`.
/// @tparam H A type satisfying the `labelled` concept.
template <typename H> requires ( requires (const H& h) { h.label ; } || requires (const H& h) { label(h) ; } )
auto getLabel(const H& h)
{
  if constexpr (requires (const H& h) { h.label ; }) return h.label;
  else return label(h);
}


/// @brief Compile-time array `{0, 1, …, N-1}` of `size_t`.
/// Used to construct `retainedAxes` covering all axes of a rank-N array, e.g. in `Master` row iteration and `superoperatorFromJump`.
template <size_t N> constexpr auto compileTimeOrdinals = [] {std::array<size_t,N> res{}; std::iota(res.begin(),res.end(),0); return res;} ();


/// @brief The multilambda pattern: combines multiple lambdas into one overload set via move construction.
/// @see https://www.scs.stanford.edu/~dm/blog/param-pack.html#multilambda
template<typename ...L>
struct multilambda : L... {
  using L::operator()...;
  constexpr multilambda(L...lambda) : L(std::move(lambda))... {}
};

} // cppqedutils


namespace boost::serialization {

/// @brief Boost.Serialization support for `cppqedutils::LogTree` (= `boost::json::object`).
/// Serializes by round-tripping through a JSON string.
/// @note Transitional: will be superseded when metadata moves to HDF5 attributes.
template <class Archive>
void serialize(Archive & ar, ::cppqedutils::LogTree& o, unsigned int version)
{
  split_free(ar, o, version);
}


template<class Archive>
void save(Archive& ar, const ::cppqedutils::LogTree& l, unsigned int) {
  ar & serialize(l) ;
}

template<class Archive>
void load(Archive& ar, ::cppqedutils::LogTree& l, unsigned int) {
  std::string s;
  ar & s ;
  l=::cppqedutils::json::parse(s).as_object();
}

} // boost::serialization