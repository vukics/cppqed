// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ElementAveraged.h"

#include "FormDouble.h"


using namespace std;



std::ostream& structure::details::streamCommon(const Averages& averages, std::ostream& os, int precision)
{
  using namespace formdouble;

  os<<'\t';
  {
    const FormDouble fd(precision);
    ranges::for_each(averages,[&](double v){os<<fd(v);});
  }
  return os;
}

