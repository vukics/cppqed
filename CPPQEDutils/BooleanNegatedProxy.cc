// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BooleanNegatedProxy.h"

#include <iostream>

/*
std::ostream& operator<<(std::ostream& os, const cpputils::BooleanNegatedProxy& b)
{
  return os<<bool(b);
}
*/

std::istream& cpputils::operator>>(std::istream& is,       cpputils::BooleanNegatedProxy& b)
{
  bool neg; is>>neg; b=neg; return is;
}
