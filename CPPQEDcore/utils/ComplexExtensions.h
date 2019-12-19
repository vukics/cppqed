// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Additional helpers for dcomp}
#ifndef CPPQEDCORE_UTILS_COMPLEXEXTENSIONS_H_INCLUDED
#define CPPQEDCORE_UTILS_COMPLEXEXTENSIONS_H_INCLUDED

#include "core_config.h"

#include <complex>

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/serialization/complex.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION


/// Double-precision complex number
/** Even though it is a type, we name it this way because we would like it to closely resemble built-in types */
typedef std::complex<double> dcomp;

/// Imaginary unit
const dcomp DCOMP_I(0,1);

inline bool isNonZero(const dcomp& c) {return bool(real(c)) || bool(imag(c));}

inline bool hasRealPart(const dcomp& c) {return bool(real(c));}
inline bool hasImagPart(const dcomp& c) {return bool(imag(c));}

inline bool  absCompare(const dcomp& c1, const dcomp& c2) {return  abs(c1)< abs(c2);}
inline bool realCompare(const dcomp& c1, const dcomp& c2) {return real(c1)<real(c2);}

#endif // CPPQEDCORE_UTILS_COMPLEXEXTENSIONS_H_INCLUDED
