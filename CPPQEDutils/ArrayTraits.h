// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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
bool isStorageContiguous(const A& a);


template<typename A>
auto size(const A& a) {return a.size();}

template<typename A>
Extents_t<Rank_v<A>> extents(const A& a);



template<typename A>
const double* data(const A& a);

template<typename A>
inline double* data(A& a) {return const_cast<double*>(data(static_cast<const A&>(a)));}



template<typename A>
A create(double* y, Extents_t<Rank_v<A>> e); ///< Create a non-owning array of data `y` with memory layout specified by `e`

template<typename A>
inline const A create(const double* y, Extents_t<Rank_v<A>> e) {return create<A>(const_cast<double*>(y),e);}

template<typename A>
A create(Extents_t<Rank_v<A>> e); ///< Owning empty array with memory layout specified by `e`

//@}



/// \name Array traversal traits
//@{

/// subscription of `a` (which might be a multi-array) with a *single* integer
template<typename A>
const auto& subscript(const A& a, size_t i) {return a[i];}


/// non-const subscription
template<typename A>
inline auto& subscript(A& a, size_t i) {return const_cast<ElementType_t<A>&>(subscript(static_cast<const A&>(a),i));}


template<typename A>
size_t subscriptLimit(const A& a) {return a.size();}

//@}


/** Specializations for std::array */

template<typename T, std::size_t N>
constexpr auto Rank_v<std::array<T,N>> = 1;


} // cppqedutils


#endif // CPPQEDCORE_UTILS_ARRAYTRAITS_H_INCLUDED
