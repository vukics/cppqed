// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BooleanNegatedProxy.h"

#include <iostream>

/*
std::ostream& operator<<(std::ostream& os, const cppqedutils::BooleanNegatedProxy& b)
{
  return os<<bool(b);
}
*/

std::istream& cppqedutils::operator>>(std::istream& is,       cppqedutils::BooleanNegatedProxy& b)
{
  bool neg; is>>neg; b=neg; return is;
}
