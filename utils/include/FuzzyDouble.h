// -*- C++ -*-
#ifndef UTILS_INCLUDE_FUZZYDOUBLE_H_INCLUDED
#define UTILS_INCLUDE_FUZZYDOUBLE_H_INCLUDED

#include "FuzzyDoubleFwd.h"

#include<boost/operators.hpp>

namespace cpputils {


class FuzzyDouble 
  : private boost::totally_ordered<FuzzyDouble>, private boost::totally_ordered<FuzzyDouble,double>
{
public:
  FuzzyDouble()         : d_( ) {}
  FuzzyDouble(double d) : d_(d) {}

  FuzzyDouble& operator=(double d) {d_=d; return *this;}

  // double& operator()() const {return d_;}

  operator double() const {return d_;}

private:
  static const double eps;

  double d_;

  friend bool operator==(FuzzyDouble, FuzzyDouble);
  friend bool operator< (FuzzyDouble, FuzzyDouble);

  friend bool operator==(FuzzyDouble, double);
  friend bool operator< (FuzzyDouble, double);
  friend bool operator> (FuzzyDouble, double);

};


} // cpputils

#endif // UTILS_INCLUDE_FUZZYDOUBLE_H_INCLUDED
