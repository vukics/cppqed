// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/** \file
 * \brief Defines wrapper functions for mathematical functions taken from libraries, and several other mathematical functions.
 * \details The principal aim of this arrangement is to localize dependence on GSL. In particular, GSL headers may be included only in .cc files.
 */
#ifndef CPPQEDCORE_UTILS_MATHEXTENSIONS_H_INCLUDED
#define CPPQEDCORE_UTILS_MATHEXTENSIONS_H_INCLUDED

#include "MathExtensionsFwd.h"

#include "ComplexExtensions.h"
#include "Exception.h"


/// Comprises wrapper functions for mathematical functions taken from libraries (Boost.Math, GSL), and several other mathematical functions.
namespace mathutils {

struct FactOverflow  : public cpputils::Exception {};
  
extern const double PI    ;
extern const double SQRTPI;
extern const double EULER ;
  
int sign(double);
int fcmp(double, double, double);
  
template<class T> const T sqr(T x) {return x*x;}

double sqr(double x);

double sqrAbs(const dcomp&);

double fact  (unsigned          ) throw(FactOverflow);
double choose(unsigned, unsigned);

template<typename T>
bool parity(T n) { return n & 1;}
// even --- false, odd --- true;


inline int round(double r) {
  return (r>0) ? floor(r+0.5) : ceil(r-0.5);
}

template<typename T>
char minusOneToThePowerOf(T n) {return parity(n) ? -1 : 1;}


/// Calculates \f$\alpha^n/\sqrt{n!}\f$ relying on the Stirling formula if \f$n\f$ is too large for explicit calculation of factorial
/**
 * The point of switching is determined by \refBoostConstruct{boost::math::max_factorial,math/doc/html/math_toolkit/factorials/sf_factorial.html}.
 *
 * The Stirling formula is a very good approximation already for moderate \f$n\f$ values:
 * \image html differencesInCoherentElements.svg
 *
 */
dcomp coherentElement(unsigned long n, const dcomp& alpha);

} // mathutils

#endif // CPPQEDCORE_UTILS_MATHEXTENSIONS_H_INCLUDED
