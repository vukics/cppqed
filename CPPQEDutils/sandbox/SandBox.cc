#include <iostream>

#include "MultiArrayComplex.h"

#include "Random.h"

#include <tuple>
//#include "MultiArrayComplex.h"

/*
using namespace cppqedutils;
using namespace boost::json;
*/

int main()
{
  auto v=boost::json::value_from(std::array{1,2,3});
  
  cppqedutils::Extents<3> /*std::array*/ a{1,2,3};
  std::cerr<<cppqedutils::toStringJSON(a)<<std::endl;

  cppqedutils::MultiArray<int,9> array{{2,4,6}};
  std::cerr<<cppqedutils::toStringJSON(array)<<std::endl;

}
