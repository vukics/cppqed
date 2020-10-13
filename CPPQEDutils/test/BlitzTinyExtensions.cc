// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzTinyExtensions.tcc"

#include <blitz/tinyvec2io.cc>

using namespace std;
using namespace blitzplusplus;


int main()
{
  {
    TinyVector<int,3> v1(0,1,2);
    TinyVector<int,4> v2(3,4,5,6);

    cout<<v1<<'&'<<v2<<" concatenated "<<concatenateTinies(v1,v2)<<" concatenated the other way round "<<concatenateTinies(v2,v1)<<endl;
  }
  {
    TinyVector<int,6> v(3,4,2,3,4,2);

    cout<<v<<" half cut "<<halfCutTiny(v)<<endl;
  }
  /* This should not compile
  {
    TinyVector<int,7> v(3,4,2,3,4,2,4);
    halfCutTiny(v);
  }
  */
#ifndef   NDEBUG
  try {
    TinyVector<int,6> v(3,4,2,1,4,2);

    cout<<"Trying to half cut "<<v<<"..."<<endl;
    halfCutTiny(v);
  }
  catch (HalfCutTinyException) {
    cout<<"HalfCutTinyException caught"<<endl; 
  }
#endif // NDEBUG

}
