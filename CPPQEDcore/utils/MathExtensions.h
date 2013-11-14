// -*- C++ -*-
/*
  Wrapper functions for mathematical functions taken from libraries,
  and several other mathematical functions.

  The principal aim of Math.h and the corresponding Math.cc is to
  localize dependence on GSL.

*/

#ifndef UTILS_MATHEXTENSIONS_H_INCLUDED
#define UTILS_MATHEXTENSIONS_H_INCLUDED

#include "MathExtensionsFwd.h"

#include "ComplexExtensions.h"
#include "Exception.h"


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


/// Calculates \f$\alpha^n/\sqrt{n!}\f$ relying on the Stirling formula if n is too large for explicit calculation of factorial
dcomp coherentElement(unsigned long n, const dcomp& alpha);

} // mathutils

#endif // UTILS_MATHEXTENSIONS_H_INCLUDED
