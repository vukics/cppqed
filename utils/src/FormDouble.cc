#include "FormDouble.h"

#include <cmath>
#include <limits>
#include <sstream>


using namespace std;


namespace {

/*
  It is easy to see that if there can be three digits in the exponent then the output takes at most 8 additional characters (together with the trailing space necessary to separate from the next field). Eg with precision=3 :
  |-1.23e-456 | => 11 characters
  |-0.0000123 |

  For positive numbers, the additional characters are at most 7.
*/


int widthPositive(int precision)
{
  return precision+int(log10(-log10(numeric_limits<double>::min())))+5;
}


int widthAny     (int precision)
{
  return widthPositive(precision)+1;
}


}



FormDouble::FormDouble(int precision) 
  : precision_(precision), width_(widthAny(precision))
{
}


const formdouble::Bound FormDouble::operator()(double d) const
{
  return formdouble::Bound(*this,d);
}


std::ostream& formdouble::operator<<(std::ostream& os, const formdouble::Bound& bf)
{
  ostringstream s; 
  // the intermediate stringstream is taken so that the manips don't affect the NEXT output
  s.precision(bf.f.precision_);
  s.width(bf.f.width_);
  s.setf(ios_base::left,ios_base::adjustfield);
  s<<bf.val;
  return os<<s.str();
}



const FormDouble formdouble::positive(int precision)
{
  const int precisionActual=max(6,precision);
  return FormDouble(precisionActual,widthPositive(precisionActual));
}


