// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Conversions.h"

#include<iostream>


using namespace std;

int main()
{
  size_t zero=0;
  size_t big=(zero-1)/10;

  cout<<"ZeroSize="<<zero<<"\nBigSize=(ZeroSize-1)/10="<<big<<endl;

  int bigInt(big);

  cout<<"int(BigSize)="<<bigInt<<endl;

  try {
    bigInt=size2Int(big);
    cout<<"int(BigSize)="<<bigInt<<endl;
  }
  catch (const boost::numeric::bad_numeric_cast& exception) {
    cout<<"Exception caught "<<exception.what()<<endl;
  }

  int negative=-1;

  try {
    big=int2Size(negative);
  }
  catch (const boost::numeric::bad_numeric_cast& exception) {
    cout<<"Exception caught "<<exception.what()<<endl;
  }

}

