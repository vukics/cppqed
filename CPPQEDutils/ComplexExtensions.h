// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Additional helpers for dcomp}
#ifndef CPPQEDCORE_UTILS_COMPLEXEXTENSIONS_H_INCLUDED
#define CPPQEDCORE_UTILS_COMPLEXEXTENSIONS_H_INCLUDED

#include <complex>

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/complex.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION


/// Double-precision complex number
/** Even though it is a type, we name it this way because we would like it to closely resemble built-in types */
typedef std::complex<double> dcomp;

using std::complex_literals::operator""i;

inline bool isNonZero(dcomp c) {return bool(real(c)) || bool(imag(c));}

inline bool hasRealPart(dcomp c) {return bool(real(c));}
inline bool hasImagPart(dcomp c) {return bool(imag(c));}

inline bool  absCompare(dcomp c1, dcomp c2) {return  abs(c1)< abs(c2);}
inline bool realCompare(dcomp c1, dcomp c2) {return real(c1)<real(c2);}

#endif // CPPQEDCORE_UTILS_COMPLEXEXTENSIONS_H_INCLUDED
