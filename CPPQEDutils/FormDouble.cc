// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "FormDouble.h"

#include <cmath>
#include <limits>


using namespace std;


/*

  It is easy to see that the output width is at most (maximal number of exponential digits) + precision + 5 characters (together with the trailing space necessary to separate from the next field).
  Eg with (maximal number of exponential digits)=3, precision=3 :
  |-1.23e-456 | => 11 characters
  |-0.0000123 |

  For positive-only numbers, the maximal width is (maximal number of exponential digits) + precision + 4

  The (maximal number of exponential digits) is expressed as
  int(log10(-log10(numeric_limits<double>::min()))) + 1

*/


int formdouble::widthPositive(int precision)
{
  return ( int(log10(-log10(numeric_limits<double>::min()))) + 1 ) + precision + 4 ;
}


int formdouble::widthAny     (int precision)
{
  return widthPositive(precision)+1;
}




FormDouble::FormDouble(int precision) 
  : FormDouble(precision,formdouble::widthAny(precision))
{
}



