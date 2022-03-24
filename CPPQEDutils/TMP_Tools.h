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


/// Cf. \refBoost{this page,fusion/doc/html/fusion/container/vector.html}
#ifndef FUSION_MAX_VECTOR_SIZE
#define FUSION_MAX_VECTOR_SIZE 20
#endif // FUSION_MAX_VECTOR_SIZE


/// Cf. \refBoost{this page,fusion/doc/html/fusion/container/list.html}
#ifndef FUSION_MAX_LIST_SIZE
#define FUSION_MAX_LIST_SIZE FUSION_MAX_VECTOR_SIZE
#endif // FUSION_MAX_LIST_SIZE

#include <boost/mpl/max_element.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/push_back.hpp>

#include <type_traits>

/// Template metaprogramming tools
namespace tmptools {

template<auto I>
using integral_c=std::integral_constant<decltype(I),I>;
  

template<typename T, bool ADD_CONST>
using ConditionalAddConst_t = std::conditional_t<ADD_CONST,std::add_const_t<T>,T>;


/// An integer \refBoostConstruct{range_c,mpl/doc/refmanual/range-c.html} starting with `Nbeg` and having `N` elements (`Nbeg ... Nbeg+N-1`)
template<unsigned N, int Nbeg> struct Range : boost::mpl::range_c<int,Nbeg,Nbeg+N> {};


/// Sequence of ordinals `0 ... N-1` based on Range
template<int N> struct Ordinals : Range<N,0> {};


/// Calculates the power \f$N_1^{N_2}\f$ @ compile time
template<unsigned N1, unsigned N2>
constexpr int Power_v=N1*Power_v<N1,N2-1>;

/** \cond SPECIALIZATION */
template<unsigned N1>
constexpr int Power_v<N1,0> =1;
/** \endcond */


/// Determines whether a compile-time sequence “numerically contains” a given value
/**
 * \tparam Seq the sequence (presumably containing integers)
 * \tparam ICW \refBoost{an integral constant wrapper,mpl/doc/refmanual/integral-constant.html}
 * 
 * The motivation for this metafunction is that Boost.MPL’s \refBoostConstruct{contains,mpl/doc/refmanual/contains.html} metafunction
 * looks for the identity of the actual types, and not the numerical equality of the contained values.
 * 
 * \see numerical_equal
 */
template<typename Seq, typename ICW>
struct numerical_contains : boost::mpl::not_<boost::is_same<typename boost::mpl::find_if<Seq,
                                                                                         boost::mpl::equal_to<boost::mpl::_,ICW> 
                                                                                         >::type,
                                                            typename boost::mpl::end<Seq>::type
                                                            >
                                             > {};

/// The `_c` version of numerical_contains, which expects a value instead of an integral constant wrapper
template<typename Seq, auto VALUE>
struct numerical_contains_c : numerical_contains<Seq,boost::mpl::integral_c<decltype(VALUE),VALUE> > {};


namespace details {

template<typename T1, typename T2>
struct value_equal : boost::mpl::bool_<T1::value==T2::value>
{};

} // details


/// Determines the numerical equality of two compile-time sequences
/**
 * \see numerical_contains for the motivation
 */
template<typename Seq1, typename Seq2>
struct numerical_equal : boost::mpl::equal<Seq1,Seq2,details::value_equal<boost::mpl::_1,boost::mpl::_2> >
{};


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

  static const auto first =N1;
  static const auto second=N2;
  
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


////////////////////////////////////////////
//
// A nonnegative compile-time integer vector
//
////////////////////////////////////////////


namespace details {

template<int V>
struct ArgumentDispatcher : boost::mpl::integral_c<int,V>
{
  static_assert( V>=0 , "Negative element in nonnegative vector" );
};

} // details

/// A non-negative compile-time vector
/**
 * It is based on and has the same characteristics as Boost.MPL’s \refBoostConstruct{vector,mpl/doc/refmanual/vector.html}.
 * 
 * Besides the check for the non-negativity of the elements @ compile time, its main motivation is that syntactic sugar that one can simply write
 * `tmptools::Vector<2,4,1,...>` instead of `boost::mpl::vector_c<int,2,4,1,...>`, where the superfluity of `int` creates considerable mental friction
 * 
 * \tparam V variadic number of integers
 * 
 */
template<int... V> 
struct Vector : boost::mpl::vector_c<int,details::ArgumentDispatcher<V>::value...> {};


typedef Vector<> V_Empty;



template<int RANK, typename V>
struct ExtendVector : boost::mpl::fold<V,
                                       V,
                                       boost::mpl::push_back<boost::mpl::_1,
                                                             boost::mpl::plus<boost::mpl::_2,
                                                                              boost::mpl::int_<RANK>
                                                                             >
                                                            >
                                      >
{};


template<int RANK, typename V>
using ExtendVector_t = typename ExtendVector<RANK,V>::type;


template <int RANK, typename V>
requires ( boost::mpl::deref<boost::mpl::max_element<V> >::type::type::value < RANK )
using NegatedView=boost::mpl::filter_view<tmptools::Ordinals<RANK>,boost::mpl::not_<tmptools::numerical_contains<V,boost::mpl::_> > >;

template <int RANK, typename V>
using NegatedVector=typename boost::mpl::copy<NegatedView<RANK,V>,boost::mpl::back_inserter<V_Empty> >::type;
// copy makes a vector of filter_view, which is important since the latter cannot be applied in all algorithms


} // tmptools



#endif // CPPQEDCORE_UTILS_TMP_TOOLS_H_INCLUDED
