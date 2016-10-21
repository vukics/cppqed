// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_FUZZYDOUBLE_H_INCLUDED
#define CPPQEDCORE_UTILS_FUZZYDOUBLE_H_INCLUDED

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

#endif // CPPQEDCORE_UTILS_FUZZYDOUBLE_H_INCLUDED
