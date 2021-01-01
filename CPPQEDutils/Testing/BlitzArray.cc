// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArray.h"

#include "Random.h"

#define BOOST_TEST_MODULE BlitzArray test
#include <boost/test/unit_test.hpp>

using namespace std;

typedef blitz::Array<dcomp,6> Array;

BOOST_AUTO_TEST_CASE( BlitzArrayTest )
{
  Array array(blitz::shape(2,4,3,2,5,4));

  cpputils::fillWithRandom<std::uniform_real_distribution<double>,std::mt19937>(array,1001);

  array.transposeSelf(1,3,0,4,2,5);

  Array arrayv(array.copy()), arrayvv(array.shape()); arrayvv=array; 
  
  cout<<array.ordering()<<arrayv.ordering()<<arrayvv.ordering()<<endl;

  dcomp* data=new dcomp[array.size()];

  memcpy(data,array.dataFirst(),array.size());

  BOOST_CHECK( all(arrayvv==array) );

}

