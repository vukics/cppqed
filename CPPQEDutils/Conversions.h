// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Some converters along the lines of \refBoost{Boost.NumericConversion,numeric/conversion/doc/html/index.html}} \todo Update & see why it is not much used in the framework
#ifndef   CPPQEDCORE_UTILS_CONVERSIONS_H_INCLUDED
#define   CPPQEDCORE_UTILS_CONVERSIONS_H_INCLUDED


#include <boost/numeric/conversion/converter.hpp>



namespace cppqedutils {


namespace conversions {


template<class T>
struct DummyRangeCheckerPolicy
{
  typedef typename T::argument_type argument_type;
  static boost::numeric::range_check_result out_of_range(argument_type) {return boost::numeric::cInRange;}
  static void validate_range(argument_type) {}
};


using OverflowHandler =
#ifndef   NDEBUG
boost::numeric::   def_overflow_handler
#else
boost::numeric::silent_overflow_handler
#endif // NDEBUG
;


template<class T>
using RangeCheckerMF = 
#ifndef   NDEBUG
  boost::numeric::UseInternalRangeChecker
#else
  DummyRangeCheckerPolicy<T>
#endif // NDEBUG
;


} // conversions



template<typename T, typename S>
struct Converter : boost::numeric::converter<T,S,
                                             boost::numeric::conversion_traits<T,S>,
                                             cppqedutils::conversions::OverflowHandler,
                                             boost::numeric::Trunc<typename boost::numeric::conversion_traits<T,S>::source_type>,
                                             boost::numeric::raw_converter<boost::numeric::conversion_traits<T,S> >,
                                             cppqedutils::conversions::RangeCheckerMF<boost::numeric::conversion_traits<T,S> >
                                             >
{};


} // cppqedutils


// from double

const cppqedutils::Converter<int,double> double2Int{};
const cppqedutils::Converter<long,double> double2Long{};
const cppqedutils::Converter<unsigned,double> double2Unsigned{};
const cppqedutils::Converter<size_t,double> double2Size_t{};
const cppqedutils::Converter<ptrdiff_t,double> double2Ptrdiff_t{};

// from long double

const cppqedutils::Converter<long,long double> longDouble2Long{};
const cppqedutils::Converter<size_t,long double> longDouble2Size_t{};
const cppqedutils::Converter<ptrdiff_t,long double> longDouble2Ptrdiff_t{};

// from size_t

const cppqedutils::Converter<double,size_t> size_t2Double{};

// from long

const cppqedutils::Converter<int,long> long2Int{};


#endif // CPPQEDCORE_UTILS_CONVERSIONS_H_INCLUDED
