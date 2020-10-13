// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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

#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>


#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/push_back.hpp>


#define DEFINE_TYPED_STATIC_CONST(typeDescription,typeName,variableName) typedef typeDescription typeName; static const typeName variableName;

#define DEFINE_INITIALIZE_TYPED_STATIC_CONST(typeDescription,typeName,variableName) typedef typeDescription typeName; static const typeName variableName=typeName();


/// Template metaprogramming tools
namespace tmptools {


/// Combines \refBoostConstruct{remove_const,type_traits/doc/html/boost_typetraits/reference/remove_const.html} and \refBoostConstruct{remove_reference,type_traits/doc/html/boost_typetraits/reference/remove_reference.html}
template<typename T>
struct RemoveConstReference : boost::remove_const<typename boost::remove_reference<T>::type>
{};


/// Applies \refBoostConstruct{add_const,type_traits/doc/html/boost_typetraits/reference/add_const.html} if `ADD_CONST = true`.
template<typename T, bool ADD_CONST>
struct ConditionalAddConst 
  : boost::mpl::eval_if_c<ADD_CONST,boost::add_const<T>,boost::mpl::identity<T> > {};


/// An integer \refBoostConstruct{range_c,mpl/doc/refmanual/range-c.html} starting with `Nbeg` and having `N` elements (`Nbeg ... Nbeg+N-1`)
template<unsigned N, int Nbeg> struct Range : boost::mpl::range_c<int,Nbeg,Nbeg+N> {};


/// Sequence of ordinals `0 ... N-1` based on Range
template<int N> struct Ordinals : Range<N,0> {};


/// Calculates the power \f$N_1^{N_2}\f$ @ compile time
template<unsigned N1, unsigned N2> struct Power       : boost::mpl::int_<N1*Power<N1,N2-1>::value> {};

/** \cond SPECIALIZATION */
template<unsigned N1>              struct Power<N1,0> : boost::mpl::int_<1>                        {};
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
template<typename Seq, typename T, T VALUE>
struct numerical_contains_c : numerical_contains<Seq,boost::mpl::integral_c<T,VALUE> > {};


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


template<bool COND, int TRUE_VALUE, int FALSE_VALUE> using integral_if_c = boost::mpl::if_c<COND,boost::mpl::int_<TRUE_VALUE>,boost::mpl::int_<FALSE_VALUE> >;

template<typename COND, int TRUE_VALUE, int FALSE_VALUE> using integral_if = integral_if_c<COND::value,TRUE_VALUE,FALSE_VALUE>;


/// Provokes a compile-time error if `N` is not even, otherwise acts as a Boost.MPL \refBoostConstruct{int_,mpl/doc/refmanual/int.html} integral constant wrapper of `N/2`
template<int N>
struct AssertEvenAndDivideBy2 : boost::mpl::int_<N/2>
{
  static_assert( 2*(N/2)==N , "Argument not even" );
};


/// A compile-time pair of integers
/** \tparam IS_EXCLUSIVE governs the exclusivity of the pair, that is, in the case of `IS_EXCLUSIVE=true`, the class provokes a compile-time error if `N1=N2` */
template<int N1, int N2, bool IS_EXCLUSIVE=true>
struct pair_c;

/// A non-exclusive pair_c that allows the members to be equal
template<int N1, int N2>
struct pair_c<N1,N2,false>
{
  //  enum { first=N1, second=N2 };

  static const int first =N1;
  static const int second=N2;
  
  /** \cond FORTESTING */
  template<int MIN, int MAX>
  struct SanityCheck
  {
    static_assert( N1>=MIN && N2>=MIN && N1<=MAX && N2<=MAX , "pair_c sanity check failed" );
  };
  /** \endcond */

};

// Apparently, in the above, the static const is DECLARED only. Some compilers (including gcc in some circumstances) need to have them DEFINED as well, which is done below. Cf EffC++ 3rd Edition, item 2.

template<int N1, int N2>
const int pair_c<N1,N2,false>::first;

template<int N1, int N2>
const int pair_c<N1,N2,false>::second;


/// An exclusive pair_c that provokes a compile-time error if the members are equal
template<int N1, int N2>
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



} // tmptools



#endif // CPPQEDCORE_UTILS_TMP_TOOLS_H_INCLUDED
