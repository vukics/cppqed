#include "BlitzTinyExtensions.h"

void f()
{
  blitz::TinyVector<int,7> v(3,4,2,3,4,2,4);
  std::cout<<blitzplusplus::halfCutTiny(v)<<std::endl;
}
