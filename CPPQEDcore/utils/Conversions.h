// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Some converters along the lines of \refBoost{Boost.NumericConversion,numeric/conversion/doc/html/index.html}} \todo Update & see why it is not much used in the framework
#ifndef   CPPQEDCORE_UTILS_CONVERSIONS_H_INCLUDED
#define   CPPQEDCORE_UTILS_CONVERSIONS_H_INCLUDED

#include "ConversionsFwd.h"

#include <boost/mpl/identity.hpp>

#include <boost/numeric/conversion/converter.hpp>



namespace cpputils {


namespace details {


template<class T>
struct DummyRangeCheckerPolicy
{
  typedef typename T::argument_type argument_type;
  static boost::numeric::range_check_result out_of_range(argument_type) {return boost::numeric::cInRange;}
  static void validate_range(argument_type) {}
};


typedef
#ifndef   NDEBUG
boost::numeric::   def_overflow_handler
#else
boost::numeric::silent_overflow_handler
#endif // NDEBUG
OverflowHandler;


template<class T> struct RangeCheckerMF : boost::mpl::identity<
#ifndef   NDEBUG
  boost::numeric::UseInternalRangeChecker
#else
  DummyRangeCheckerPolicy<T>
#endif // NDEBUG
>
{};


} // details



template<typename T, typename S>
struct Converter : boost::numeric::converter<T,S,
                                             boost::numeric::conversion_traits<T,S>,
                                             cpputils::details::OverflowHandler,
                                             boost::numeric::Trunc<typename boost::numeric::conversion_traits<T,S>::source_type>,
                                             boost::numeric::raw_converter<boost::numeric::conversion_traits<T,S> >,
                                             typename cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<T,S> >::type
                                             >
{};


} // cpputils


// Double2Int

typedef boost::numeric::converter<int,
                                  double,
                                  boost::numeric::conversion_traits<int,double>,
                                  cpputils::details::OverflowHandler,
                                  boost::numeric::Trunc<boost::numeric::conversion_traits<int,double>::source_type>,
                                  boost::numeric::raw_converter<boost::numeric::conversion_traits<int,double> >,
                                  cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<int,double> >::type
                                  > Double2Int;

const Double2Int double2Int=Double2Int();


// Long2Int

typedef boost::numeric::converter<int,
                                  long,
                                  boost::numeric::conversion_traits<int,long>,
                                  cpputils::details::OverflowHandler,
                                  boost::numeric::Trunc<boost::numeric::conversion_traits<int,long>::source_type>,
                                  boost::numeric::raw_converter<boost::numeric::conversion_traits<int,long> >,
                                  cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<int,long> >::type
                                  > Long2Int;

const Long2Int long2Int=Long2Int();


// Size2Int

typedef boost::numeric::converter<int,
                                  size_t,
                                  boost::numeric::conversion_traits<int,size_t>,
                                  cpputils::details::OverflowHandler,
                                  boost::numeric::Trunc<boost::numeric::conversion_traits<int,size_t>::source_type>,
                                  boost::numeric::raw_converter<boost::numeric::conversion_traits<int,size_t> >,
                                  cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<int,size_t> >::type
                                  > Size2Int;

const Size2Int size2Int=Size2Int();


// Idx2Int

typedef boost::numeric::converter<int,
                                  ptrdiff_t,
                                  boost::numeric::conversion_traits<int,ptrdiff_t>,
                                  cpputils::details::OverflowHandler,
                                  boost::numeric::Trunc<boost::numeric::conversion_traits<int,ptrdiff_t>::source_type>,
                                  boost::numeric::raw_converter<boost::numeric::conversion_traits<int,ptrdiff_t> >,
                                  cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<int,ptrdiff_t> >::type
                                  > Idx2Int;

const Idx2Int idx2Int=Idx2Int();


// Int2Size

typedef boost::numeric::converter<size_t,
                                  int,
                                  boost::numeric::conversion_traits<size_t,int>,
                                  cpputils::details::OverflowHandler,
                                  boost::numeric::Trunc<boost::numeric::conversion_traits<size_t,int>::source_type>,
                                  boost::numeric::raw_converter<boost::numeric::conversion_traits<size_t,int> >,
                                  cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<size_t,int> >::type
                                  > Int2Size;

const Int2Size int2Size=Int2Size();


// Size2Double

typedef boost::numeric::converter<double,
                                  size_t,
                                  boost::numeric::conversion_traits<double,size_t>,
                                  cpputils::details::OverflowHandler,
                                  boost::numeric::Trunc<boost::numeric::conversion_traits<double,size_t>::source_type>,
                                  boost::numeric::raw_converter<boost::numeric::conversion_traits<double,size_t> >,
                                  cpputils::details::RangeCheckerMF<boost::numeric::conversion_traits<double,size_t> >::type
                                  > Size2Double;

const Size2Double size2Double=Size2Double();

#endif // CPPQEDCORE_UTILS_CONVERSIONS_H_INCLUDED
