// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Combinatorics.h"

#include "MathExtensions.h"


using namespace std;

using blitz::Range;


cppqedutils::CWR_Dir::CWR_Dir(size_t n, size_t k)
  : impl_(recurse(Impl(choose(n+k-1,k),n),k))
{}


const cppqedutils::CWR_Dir::Configuration cppqedutils::CWR_Dir::operator[](size_t i) const
{
  return impl_(i,Range::all());
}



const cppqedutils::CWR_Dir::Impl cppqedutils::CWR_Dir::recurse(cppqedutils::CWR_Dir::Impl dir, size_t k)
{
  const size_t 
    n=dir.extent(1); // Current n

  if (n==1 || !k) dir=k;
  else 
    for (int
           putHere=k,
           extentHere=1,
           topLeft=0;
         putHere>=0;
         (putHere--,
          topLeft+=extentHere,
          extentHere=choose((n-1)+(k-putHere)-1,(k-putHere))
          )) {
      
      Range zeroth(topLeft,topLeft+extentHere-1);
      dir(zeroth,0)=putHere;
      recurse(dir(zeroth,Range(1,n-1)),k-putHere);
    }

  return dir;

}

