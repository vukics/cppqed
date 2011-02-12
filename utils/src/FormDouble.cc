#include "FormDouble.h"

#include<limits>
#include<sstream>



using namespace std;


namespace formdouble {


FormDouble::FormDouble(int precision) 
  : precision_(precision), width_(precision+(std::numeric_limits<double>::min()<1e-99 ? 8 : 7))
    /*
      It is easy to see that if there can be three digits in the exponent
      then the output takes at most 8 additional characters.
      (together with the trailing space necessary to separate from the next field)
      Eg with precision=3
      |-1.23e-456 | => 11 characters
      |-0.0000123 |
    */
{
}


const BoundFormDouble FormDouble::operator()(double d) const
{
  return BoundFormDouble(*this,d);
}


std::ostream& operator<<(std::ostream& os, const BoundFormDouble& bf)
{
  ostringstream s; 
  // the intermediate stringstream is taken so that the manips don't affect the NEXT output
  s.precision(bf.f.precision_);
  s.width(bf.f.width_);
  s.setf(ios_base::left,ios_base::adjustfield);
  s<<bf.val;
  return os<<s.str();
}


} // formdouble

