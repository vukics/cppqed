// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArray.h"

#include "Randomized.h"

#include "Range.h"

#define BOOST_TEST_MODULE BlitzArray test
#include <boost/test/unit_test.hpp>

#include<boost/bind.hpp>

using namespace std;
using namespace randomized;

typedef blitz::Array<dcomp,6> Array;

typedef blitzplusplus::Array<dcomp,6> MyArray;


BOOST_AUTO_TEST_CASE( BlitzArrayTest )
{
  Array array(blitz::shape(2,4,3,2,5,4));

  fillWithRandom(array);

  array.transposeSelf(1,3,0,4,2,5);

  Array arrayv(array.copy()), arrayvv(array.shape()); arrayvv=array; 
  
  cout<<array.ordering()<<arrayv.ordering()<<arrayvv.ordering()<<endl;

  dcomp* data=new dcomp[array.size()];

  memcpy(data,array.dataFirst(),array.size());

  blitzplusplus::Array<dcomp,6> myArray(MyArray(array,blitzplusplus::ShallowCopy()).clone(data),blitzplusplus::ShallowCopy());

  BOOST_CHECK( all(arrayvv==array) && all((2.*myArray)!=array) );

  // return 0;

}

