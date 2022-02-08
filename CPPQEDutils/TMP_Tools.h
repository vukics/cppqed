// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Template metaprogramming tools, extending (and based on) Boost.MPL.}
#ifndef CPPQEDCORE_UTILS_TMP_TOOLS_H_INCLUDED
#define CPPQEDCORE_UTILS_TMP_TOOLS_H_INCLUDED

/// The largest-rank `blitz::Array` for which the mixed-mode subscripting can be used.
/**
 * A value larger than 11 will work only if our own Blitz++ version is used, where the indexing member functions and constructors are
 * generated through preprocessor metaprogramming.
 *
 * For our version of Blitz++, cf. http://sourceforge.net/p/cppqed/blitz/ci/default/tree/
 * 
 */
#ifndef BLITZ_ARRAY_LARGEST_RANK
#define BLITZ_ARRAY_LARGEST_RANK 11
#endif // BLITZ_ARRAY_LARGEST_RANK


#include <type_traits>

#include <boost/hana.hpp>

namespace hana=boost::hana;
#include <boost/mpl/deref.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/filter_view.hpp>


/// Template metaprogramming tools
namespace tmptools {

template<auto I>
using integral_c=std::integral_constant<decltype(I),I>;
  

template<typename T, bool ADD_CONST>
using ConditionalAddConst_t = std::conditional_t<ADD_CONST,std::add_const_t<T>,T>;


/// Calculates the power \f$N_1^{N_2}\f$ @ compile time
template<unsigned N1, unsigned N2>
constexpr int Power_v=N1*Power_v<N1,N2-1>;

template<unsigned N1>
constexpr int Power_v<N1,0> =1;


template<bool COND, auto TRUE_VALUE, auto FALSE_VALUE>
constexpr auto integral_if_c_v = std::conditional_t<COND,integral_c<TRUE_VALUE>,integral_c<FALSE_VALUE>>::value;

template<typename COND, auto TRUE_VALUE, auto FALSE_VALUE>
constexpr auto integral_if_v = integral_if_c_v<COND::value,TRUE_VALUE,FALSE_VALUE>;


/// A compile-time pair of integral types
/** \tparam IS_EXCLUSIVE governs the exclusivity of the pair, that is, in the case of `IS_EXCLUSIVE=true`, the class provokes a compile-time error if `N1=N2` */
template<auto N1, auto N2, bool IS_EXCLUSIVE=true>
struct pair_c;

/// A non-exclusive pair_c that allows the members to be equal
template<auto N1, auto N2>
struct pair_c<N1,N2,false>
{
  //  enum { first=N1, second=N2 };

  static constexpr auto first =N1;
  static constexpr auto second=N2;
  
  /** \cond FORTESTING */
  template<auto MIN, auto MAX>
  struct SanityCheck
  {
    static_assert( N1>=MIN && N2>=MIN && N1<=MAX && N2<=MAX , "pair_c sanity check failed" );
  };
  /** \endcond */

};

// // Apparently, in the above, the static const is DECLARED only. Some compilers (including gcc in some circumstances) need to have them DEFINED as well, which is done below. Cf EffC++ 3rd Edition, item 2.
// 
// template<int N1, int N2>
// const int pair_c<N1,N2,false>::first;
// 
// template<int N1, int N2>
// const int pair_c<N1,N2,false>::second;


/// An exclusive pair_c that provokes a compile-time error if the members are equal
template<auto N1, auto N2>
struct pair_c<N1,N2,true> : pair_c<N1,N2,false> 
{
  static_assert( N1!=N2 , "pair_c with equal elements" );
};


template<int... i>
constexpr auto vector = hana::tuple_c<int,i...>;


constexpr auto vEmpty = vector<>;


template<int... i>
using Vector = decltype(vector<i...>);


using V_Empty = Vector<>;


template<int begin, int end>
constexpr auto range = hana::range_c<int,begin,end>;


template<int begin, int end>
using Range = decltype(range<begin,end>);


template<int end>
constexpr auto ordinals = range<0,end>;


template<int end>
using Ordinals = decltype(ordinals<end>);



template <typename Sequence>
using CopyToVector=typename boost::mpl::copy<Sequence,boost::mpl::back_inserter<V_Empty> >::type;


template <int RANK, typename V>
requires ( boost::mpl::deref<boost::mpl::max_element<V> >::type::type::value < RANK )
using NegatedView=boost::mpl::filter_view<tmptools::Ordinals<RANK>,boost::mpl::not_<tmptools::numerical_contains<V,boost::mpl::_> > >;

template <int RANK, typename V>
using NegatedVector=CopyToVector<NegatedView<RANK,V>>;
// copy makes a vector of filter_view, which is important since the latter cannot be applied in all algorithms

} // tmptools



#endif // CPPQEDCORE_UTILS_TMP_TOOLS_H_INCLUDED
