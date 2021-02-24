// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Declarations of traits functions for adapting array types to generic functions}
#ifndef   CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED
#define   CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED

#include <array>
#include <optional>
#include <vector>

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


template<typename A>
bool isStorageContiguous(const A& a);


template<typename A>
auto size(const A& a) {return a.size();}

template<typename A>
std::vector<size_t> dimensions(const A& a);



template<typename A>
const double* data(const A& a);

template<typename A>
      double* data(      A& a);



template<typename A>
      A create(      double* y, const A& a); ///< Clone (create a non-owning array of data `y` of the same memory layout as `a`)

template<typename A>
const A create(const double* y, const A& a); ///< Const clone (create a const non-owning array of data `y` of the same memory layout as `a`)


template<typename A>
A create(const A& a); ///< Empty clone (create a newly allocated owning empty array of the same memory layout as `a`)
//@}



/// \name Array traversal traits
//@{

/// subscription of `a` (which might be a multi-array) with a *single* integer
template<typename A>
const auto& subscript(const A& a, size_t i) {return a[i];}


/// non-const subscription
template<typename A>
auto& subscript(A& a, size_t i) {return a[i];}


template<typename A>
size_t subscriptLimit(const A& a) {return a.size();}

//@}


/** Specializations for std::array */

template<typename T, std::size_t N>
constexpr auto Rank_v<std::array<T,N>> = 1;


} // cppqedutils


#endif // CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED
