#include "MultiArray.h"

#include <iostream>
#include <tuple>

int main()
{
  cppqedutils::MultiArray<double,5> ma;
  
  double& v1=ma(1,5,6,3,4) ;
  // double& v2=ma(1,5,6,std::tuple<int>{3},4) ;

  // cppqedutils::MultiArray<const double,5> mac; // "std::vector must have a non-const, non-volatile value_type"

  cppqedutils::MultiArrayView<const double,5> mavc;

  const double& v3=mavc(1,5,6,3,4) ;
  // double& v4=mavc(1,5,6,3,4) ;
  
  // double& v5=ma(1,5,6,3) ;
  // double& v6=ma(1,5,6,3,4,2) ;
  
}
