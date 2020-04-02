// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "KeyPrinter.h"

#include <boost/range/algorithm/for_each.hpp>

#include <boost/lambda/lambda.hpp>

#include <iostream>
#include <iomanip>

using namespace std;


ostream& cpputils::KeyPrinter::displayKey(ostream& os, size_t& i) const
{
  namespace bll=boost::lambda;

  os<<"# "<<keyTitle_;
  boost::for_each(keyLabels_,os<<bll::constant("\n# ")<<bll::constant(setw(2))<<bll::var(i)++<<". "<<bll::_1);
  return os<<endl;
}

