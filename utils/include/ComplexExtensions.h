// -*- C++ -*-
#ifndef UTILS_INCLUDE_COMPLEXEXTENSIONS_H_INCLUDED
#define UTILS_INCLUDE_COMPLEXEXTENSIONS_H_INCLUDED

#include<complex>

#ifndef DO_NOT_USE_BOOST_SERIALIZATION
#include <boost/serialization/complex.hpp>
#endif // DO_NOT_USE_BOOST_SERIALIZATION


typedef std::complex<double> dcomp;
// Even though it is a type we name it this way because we would like it to closely resemble built-in types

const dcomp DCOMP_I(0,1);

inline bool isNonZero(const dcomp& c) {return bool(real(c)) || bool(imag(c));}

inline bool hasRealPart(const dcomp& c) {return bool(real(c));}
inline bool hasImagPart(const dcomp& c) {return bool(imag(c));}

inline bool  absCompare(const dcomp& c1, const dcomp& c2) {return  abs(c1)< abs(c2);}
inline bool realCompare(const dcomp& c1, const dcomp& c2) {return real(c1)<real(c2);}

#endif // UTILS_INCLUDE_COMPLEXEXTENSIONS_H_INCLUDED