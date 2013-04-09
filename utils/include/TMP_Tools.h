// -*- C++ -*-
#ifndef UTILS_INCLUDE_TMP_TOOLS_H_INCLUDED
#define UTILS_INCLUDE_TMP_TOOLS_H_INCLUDED

#ifndef BLITZ_ARRAY_LARGEST_RANK
#define BLITZ_ARRAY_LARGEST_RANK 11
#endif // BLITZ_ARRAY_LARGEST_RANK

#ifndef FUSION_MAX_VECTOR_SIZE
#define FUSION_MAX_VECTOR_SIZE 20
#endif // FUSION_MAX_VECTOR_SIZE

#ifndef FUSION_MAX_LIST_SIZE
#define FUSION_MAX_LIST_SIZE FUSION_MAX_VECTOR_SIZE
#endif // FUSION_MAX_LIST_SIZE


#include "TMP_ToolsFwd.h"

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

#include <boost/mpl/assert.hpp>

#define DEFINE_TYPED_STATIC_CONST(typeDescription,typeName,variableName) typedef typeDescription typeName; static const typeName variableName;

#define DEFINE_INITIALIZE_TYPED_STATIC_CONST(typeDescription,typeName,variableName) typedef typeDescription typeName; static const typeName variableName=typeName();


namespace tmptools {


template<typename T>
struct RemoveConstReference : boost::remove_const<typename boost::remove_reference<T>::type>
{};


template<typename T, bool CONST>
struct ConditionalAddConst 
  : boost::mpl::eval_if_c<CONST,boost::add_const<T>,boost::mpl::identity<T> > {};


namespace details {

typedef boost::mpl::range_c<int,0,0> emptyRange;

} // details


template<int N, int Nbeg> struct Range : boost::mpl::range_c<int,Nbeg,Nbeg+N>
{
  BOOST_MPL_ASSERT_MSG( N >= 0, RANGE_with_NEGATIVE_LENGTH, (boost::mpl::int_<N>) );
};


template<int N> struct Ordinals : Range<N,0> {};


template<int N1, int N2> struct Power       : boost::mpl::int_<N1*Power<N1,N2-1>::value> {};

template<int N1>         struct Power<N1,0> : boost::mpl::int_<1>                        {};


template<typename Seq, typename ICW>
struct numerical_contains : boost::mpl::not_<boost::is_same<typename boost::mpl::find_if<Seq,
                                                                                         boost::mpl::equal_to<boost::mpl::_,ICW> 
                                                                                         >::type,
                                                            typename boost::mpl::end<Seq>::type
                                                            >
                                             > {};


template<typename Seq, typename T, T VALUE>
struct numerical_contains_c : numerical_contains<Seq,boost::mpl::integral_c<T,VALUE> > {};


namespace details {

template<typename T1, typename T2>
struct value_equal : boost::mpl::bool_<T1::value==T2::value>
{};

} // details


template<typename Seq1, typename Seq2>
struct numerical_equal : boost::mpl::equal<Seq1,Seq2,details::value_equal<boost::mpl::_1,boost::mpl::_2> >
{};



template<int N>
struct IsEvenAssert : boost::mpl::int_<N/2>
{
  BOOST_MPL_ASSERT_MSG( (2*(N/2)==N) , ARGUMENT_NOT_EVEN, (boost::mpl::int_<N>) );
};


template<int N1, int N2, bool IS_EXCLUSIVE=true>
struct pair_c;


template<int N1, int N2>
struct pair_c<N1,N2,false>
{
  //  enum { first=N1, second=N2 };

  static const int first =N1;
  static const int second=N2;
  
  template<int MIN, int MAX>
  struct SanityCheck
  {
    BOOST_MPL_ASSERT_MSG( (N1>=MIN && N2>=MIN && N1<=MAX && N2<=MAX) , PAIR_C_SANITY_CHECK_FAILED , (SanityCheck) );
  };

};

// Apparently, in the above, the static const is DECLARED only. Some compilers (including gcc in some circumstances) need to have them DEFINED as well, which is done below. Cf EffC++ 3rd Edition, item 2.

template<int N1, int N2>
const int pair_c<N1,N2,false>::first;

template<int N1, int N2>
const int pair_c<N1,N2,false>::second;


template<int N1, int N2>
struct pair_c<N1,N2,true> : pair_c<N1,N2,false> 
{
  BOOST_MPL_ASSERT_MSG( (N1!=N2) , PAIR_C_with_EQUAL_ELEMENTS , (pair_c) );  
};


////////////////////////////////////////////
//
// A nonnegative compile-time integer vector
//
////////////////////////////////////////////


template<int V>
struct ArgumentDispatcher : boost::mpl::integral_c<int,V>
{
  BOOST_MPL_ASSERT_MSG( V>=0 , NEGATIVE_ELEMENT_in_NONNEGATIVE_VECTOR, (boost::mpl::int_<V>) );
};


// Workaround for gcc <4.7 http://stackoverflow.com/a/11297765
template<template <typename...> class T, typename... Args>
struct Join : boost::mpl::identity<T<Args...> > 
{};


template<int... V> 
struct Vector : // boost::mpl::vector_c<int,ArgumentDispatcher<V>::value...>
  Join<boost::mpl::vector,ArgumentDispatcher<V>...>::type
{};


typedef Vector<> V_Empty;


} // tmptools



#endif // UTILS_INCLUDE_TMP_TOOLS_H_INCLUDED
