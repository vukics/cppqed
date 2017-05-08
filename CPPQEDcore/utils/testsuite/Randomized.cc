// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArray.h"
#include "Randomized.h"

#define BOOST_TEST_MODULE Randomized test
#include <boost/test/unit_test.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <fstream>

using namespace randomized;
using namespace std;


const size_t seed=1001;

BOOST_AUTO_TEST_CASE( SerializationTest )
{
  DArray<1> array(10), arrayInterrupted(10), firstHalf(arrayInterrupted(blitz::Range(0,4))), secondHalf(arrayInterrupted(blitz::Range(5,9)));
  
  fillWithRandom(array,seed);
  
  Randomized::Ptr ranFirst(MakerGSL()(seed)), ranSecond(MakerGSL()(seed));
  
  fillWithRandom(firstHalf,ranFirst);
  
  {
    ofstream ofs("/tmp/C++QED_utils_testsuite_Randomized.d");
    boost::archive::binary_oarchive ar(ofs);
    ar & *ranFirst.get();
  }

  // ... some time later restore the other Randomized object instance to its state
  {
    ifstream ifs("/tmp/C++QED_utils_testsuite_Randomized.d");
    boost::archive::binary_iarchive ar(ifs);
    ar & *ranSecond.get();
  }

  fillWithRandom(secondHalf,ranSecond);
  
  BOOST_CHECK( all(array==arrayInterrupted) );
  
}
