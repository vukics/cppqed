// -*- C++ -*-
#ifndef _TEMPLATE_METAPROG_TOOLS_H
#define _TEMPLATE_METAPROG_TOOLS_H

#include "TMP_ToolsFwd.h"

#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>


#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/equal.hpp>

#include <boost/mpl/assert.hpp>

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>


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
struct IsEvenAssert
{
  BOOST_MPL_ASSERT_MSG( (2*(N/2)==N) , ARGUMENT_NOT_EVEN, (boost::mpl::int_<N>) );
  static const int value=N/2;
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

/*
  While the adopted solution relying on a nested series of
  mpl::eval_if-s may seem formidable, it provides the syntactic sugar
  of being able to write 

  Vector<1,3,2> 

  which immediately produces the Vector with the desired elements.
  Due to the lazy evaluation of mpl::eval_if, the solution should not
  be very inefficient.

  NEEDS_WORK Check Abrahams-Gurtuvoy if there is a better solution for this.
*/

#define TMPTOOLS_VECTOR_DEFAULT_ARG 10011001

template<int V, int N>
struct ArgumentDispatcher : boost::mpl::bool_<(V==TMPTOOLS_VECTOR_DEFAULT_ARG)>
{
  BOOST_MPL_ASSERT_MSG( (V>=0 || V==TMPTOOLS_VECTOR_DEFAULT_ARG) , NEGATIVE_ELEMENT_in_NONNEGATIVE_VECTOR_at, (boost::mpl::int_<N>) );
};


#define VECTOR_MAX_SIZE 20
#define VECTOR_print(z,n,unused) typename boost::mpl::eval_if<ArgumentDispatcher<V##n,n>,boost::mpl::vector_c<int,BOOST_PP_ENUM_PARAMS(n,V)> BOOST_PP_COMMA_IF(n)
#define VECTOR_trailingPrint(z,n,unused) >::type


template<BOOST_PP_ENUM_BINARY_PARAMS(VECTOR_MAX_SIZE,int V,=TMPTOOLS_VECTOR_DEFAULT_ARG BOOST_PP_INTERCEPT)> 
struct Vector
  : boost::mpl::eval_if_c<V0==TMPTOOLS_VECTOR_DEFAULT_ARG,boost::mpl::vector_c<int>,
			  BOOST_PP_REPEAT_FROM_TO(1,VECTOR_MAX_SIZE,VECTOR_print,~)
                          boost::mpl::vector_c<int,BOOST_PP_ENUM_PARAMS(VECTOR_MAX_SIZE,V)>
                          BOOST_PP_REPEAT(VECTOR_MAX_SIZE,VECTOR_trailingPrint,~)
{};


#undef  VECTOR_trailingPrint
#undef  VECTOR_print
#undef  VECTOR_MAX_SIZE


typedef Vector<> V0;


} // tmptools



#endif // _TEMPLATE_METAPROG_TOOLS_H
