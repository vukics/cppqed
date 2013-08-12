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
    re=pow((alpha/sqrt(n/2.71828)),n);
  }

  return re;
}
