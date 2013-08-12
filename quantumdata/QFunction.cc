#include "QFunction.h"

#include <boost/math/special_functions/factorials.hpp>

using namespace mathutils;
using namespace boost::math;

double quantumdata::details::qFunctionHelper(size_t n, const dcomp& alpha)
{
  if (n<100) {
    re=pow(alpha, n)/sqrt(factorial<long double>(n));
  }
  else {
    re=(1./(sqrt(sqrt(2.*3.1415*n))))*pow((alpha/sqrt(n/exp(1.))),n);
  }

  return re;
}
