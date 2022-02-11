// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Declarations of traits functions for adapting array types to generic functions}
#ifndef   CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED
#define   CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED

#include <array>
#include <optional>
#include <stdexcept>


namespace cppqedutils {


/// template metafunction for the rank (arity) of the multi-array `A`
template<typename>
constexpr auto Rank_v=std::nullopt;

/// template metafunction providing an identifier string for the multi-array `A`
template<typename>
constexpr auto TypeID_v=std::nullopt;

/// template metafunction returning (by convention, as a member typedef `type`) the type of elements of the multi-array `A`
template<typename A>
struct ElementType
{
  using type=typename A::value_type;
};

template<typename A>
using ElementType_t = typename ElementType<A>::type;

/// \name Array memory traits
//@{


template <size_t Rank>
using Extents_t=std::array<size_t,Rank>;


struct NonContiguousStorageException : public std::invalid_argument {using std::invalid_argument::invalid_argument;};


template<typename A>
struct IsStorageContiguous;


template<typename A>
struct Size
{
  static auto _(const A& a) {return a.size();}
};


template<typename A>
struct Extents;



template<typename A>
struct Data_c;


template<typename A>
struct Data
{
  static double* _(A& a) {return const_cast<double*>(Data_c<A>::_(a));}
};


template<typename A>
struct Create;


template<typename A>
struct Create_c
{
  /// Create a non-owning array of data `y` with memory layout specified by `e`
  static const A _(const double* y, Extents_t<Rank_v<A>> e) {return Create<A>::_(const_cast<double*>(y),e);}
};


///< Owning empty array with memory layout specified by Extents
template<typename A>
struct CreateFromExtents;


template<typename A>
struct Copy
{
  /// Copy of a with its own memory
  static A _(const A& a) {return a;}
};

//@}



/// \name Array traversal traits
//@{

/// subscription of `a` (which might be a multi-array) with a *single* integer
template<typename A>
struct Subscript_c
{
  static const auto& _(const A& a, size_t i) {return a[i];}
};


/// non-const subscription
template<typename A>
struct Subscript
{
  static auto& _(A& a, size_t i) {return const_cast<ElementType_t<A>&>(Subscript_c<A>::_(a,i));}
};


template<typename A>
struct SubscriptLimit
{
  static auto _(const A& a) {return a.size();}
};


template<typename A>
struct Stride;


//@}


/** Specializations for std::array */

template<typename T, std::size_t N>
constexpr auto Rank_v<std::array<T,N>> = 1;


} // cppqedutils


#endif // CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED
