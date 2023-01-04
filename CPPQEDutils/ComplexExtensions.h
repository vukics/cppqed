// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Additional helpers for dcomp}
#pragma once

#include <complex>

#include <boost/serialization/complex.hpp>


/// Double-precision complex number
/** Even though it is a type, we name it this way because we would like it to closely resemble built-in types */
typedef std::complex<double> dcomp;

using std::complex_literals::operator""i;

inline bool isNonZero(dcomp c) {return bool(real(c)) || bool(imag(c));}

inline bool hasRealPart(dcomp c) {return bool(real(c));}
inline bool hasImagPart(dcomp c) {return bool(imag(c));}

inline bool  absCompare(dcomp c1, dcomp c2) {return  abs(c1)< abs(c2);}
inline bool realCompare(dcomp c1, dcomp c2) {return real(c1)<real(c2);}

/// Generic sqr from this Q&A: https://stackoverflow.com/a/64849248/1171157
auto sqr(auto&& x) // return type is non-reference or trailing
noexcept(noexcept(x*x)) // propagate noexcept
-> decltype(x*x) // enable SFINAE
{ return x * x; }


double relativeDeviation(const auto& a, const auto& b) {return std::abs(a-b)/(std::abs(a)+std::abs(b));}
