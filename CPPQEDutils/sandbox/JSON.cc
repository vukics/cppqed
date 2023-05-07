#include "Traits.h"

#include <iostream>


using cppqedutils::LogTree;

int main()
{
  LogTree anObject{{"aDouble",1.32},{"anInt",-31},{"aString","something"}};

  assert(anObject.is_object());

  LogTree nestedObject{{"anStdArray",std::array{3,2,5}},{"memberObject",anObject}};

  assert(nestedObject.is_object());

  LogTree anArray{anObject,nestedObject,"anotherSomethig",{"anotherDouble",3.14},3.14};

  assert(anArray.is_array());

  std::cerr<<std::setw(2)<<anObject<<std::setw(2)<<nestedObject<<std::setw(2)<<anArray<<std::endl;
}
