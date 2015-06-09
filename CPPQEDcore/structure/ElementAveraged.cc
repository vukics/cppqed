// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ElementAveraged.h"

#include "FormDouble.tcc"

#include <boost/range/algorithm/for_each.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


using namespace std;



std::ostream& structure::details::displayCommon(const AveragedCommon::Averages& averages, std::ostream& os, int precision)
{
  using namespace formdouble;

  namespace bll=boost::lambda;

  os<<'\t';
  {
    const FormDouble fd(precision);
    boost::for_each(averages,os<<bll::bind(&FormDouble::operator()<double>,&fd,bll::_1));
  }
  return os;
}

