// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ElementAveraged.h"

#include "FormDouble.h"

#include <boost/range/algorithm/for_each.hpp>


using namespace std;



std::ostream& structure::details::streamCommon(const AveragedCommon::Averages& averages, std::ostream& os, int precision)
{
  using namespace formdouble;

  os<<'\t';
  {
    const FormDouble fd(precision);
    boost::for_each(averages,[&](double v){os<<fd(v);});
  }
  return os;
}

