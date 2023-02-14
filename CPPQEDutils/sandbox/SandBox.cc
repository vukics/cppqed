#include <iostream>

#include "MultiArrayComplex.h"

#include "Random.h"

#include <tuple>
//#include "MultiArrayComplex.h"

using namespace cppqedutils;
using namespace boost::json;

int main()
{
  auto v=value_from(std::array{1,2,3});
  
  Extents<3> /*std::array*/ a{1,2,3};
  std::cerr<<toStringJSON(a)<<std::endl;
}
