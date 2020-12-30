// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Helpers for `blitz::TinyVector`s}
#ifndef CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_H_INCLUDED
#define CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_H_INCLUDED

#include "Conversions.h"
#include "TMP_Tools.h"

#include <blitz/array.h> // this is needed in order to pull in the expression template mechanism

#include <boost/mpl/for_each.hpp>

#ifndef   NDEBUG
#include <stdexcept>
#endif // NDEBUG

namespace blitzplusplus {

/// Concatenates `tiny1` and `tiny2`
/**
 * Leverages the compile-time–runtime facility \refBoostConstruct{for_each,mpl/doc/refmanual/for-each.html} from Boost.MPL.
 * 
 * \tparam T1 basic type of one of the operand tiny vectors *and the result*
 * \tparam T2 basic type of the other operand tiny vector. Must be convertible to T1.
 * 
 */
template<typename T1, typename T2, int RANK1, int RANK2>
auto
concatenateTinies(const blitz::TinyVector<T1,RANK1>& op1, const blitz::TinyVector<T2,RANK2>& op2)
{
  blitz::TinyVector<T1,RANK1+RANK2> res;

  boost::mpl::for_each<tmptools::Ordinals<RANK1> >([&](auto t) {
    static const auto i=decltype(t)::value;
    res(i)=op1(i);});
  boost::mpl::for_each<tmptools::Ordinals<RANK2> >([&](auto t) {
    static const auto i=decltype(t)::value;
    res(i+RANK1)=cpputils::Converter<T1,T2>::convert(op2(i));});

  return res;
    
}


/// Returns the first half of a tiny vector containing two equal halves
/**
 * In debug mode, the equality of the two halves is checked, and a HalfCutTinyException is thrown @ violation.
 * 
 * \tparam T basic type of the operand tiny vector
 * \tparam TWO_TIMES_RANK length of the ”
 * 
 */
template<typename T, int TWO_TIMES_RANK>
std::enable_if_t<TWO_TIMES_RANK%2==0,
                 blitz::TinyVector<T,TWO_TIMES_RANK/2>>
halfCutTiny(const blitz::TinyVector<T,TWO_TIMES_RANK>& tiny)
{
  static const int RANK=TWO_TIMES_RANK/2;

#ifndef   NDEBUG
  boost::mpl::for_each<tmptools::Ordinals<RANK> >([&](auto arg) {
    using U = decltype(arg);
    if (tiny(U::value)!=tiny(U::value+RANK)) throw std::invalid_argument("In halfCutTiny");
  });
#endif // NDEBUG

  blitz::TinyVector<T,RANK> res;

  boost::mpl::for_each<tmptools::Ordinals<RANK> >([&](auto t) {
    static const auto i=decltype(t)::value;
    res(i)=tiny(i);});

  return res;

}


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_H_INCLUDED
