#include "Traits.h"

#include <array>
#include <iostream>


using namespace cppqedutils ;

int main()
{
  LogTree anObject{{"aDouble",1.32},{"anInt",-31},{"aString","something"}};

  // assert(anObject.is_object());

  LogTree nestedObject{{"anStdArray",json::value_from(std::array{3,2,5}).as_array()},{"memberObject",anObject}};

  // assert(nestedObject.is_object());

  json::array anArray{anObject,nestedObject,"anotherSomethig",{"anotherDouble",3.14},3.14};

  // assert(anArray.is_array());

  std::cerr<<anObject<<nestedObject<<anArray<<std::endl;
}
