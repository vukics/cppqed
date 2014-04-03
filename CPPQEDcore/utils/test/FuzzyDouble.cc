// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "FuzzyDouble.h"

#include<iostream>
#include<limits>

using namespace std;
using namespace cpputils;

int main()
{
  {
  double a=3.14, b=3.141;

  FuzzyDouble aF(a), bF(b);

  // aF-=FuzzyDouble(1.);

  cout<<aF+bF<<endl;
  cout<<a +bF<<endl;

  cout<<(a==aF)<<endl;
  cout<<(b>aF)<<endl;
  cout<<(b>=aF)<<endl;
  cout<<(b<aF)<<endl;

  cout<<-aF<<endl;

  }

  cout<<(3*1e-4==0.0003)<<endl;
  cout<<(FuzzyDouble(3*1e-4)==0.0003)<<endl;

  cout<<numeric_limits<double>::epsilon()<<endl;

}
