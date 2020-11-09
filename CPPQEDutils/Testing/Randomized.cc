// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArray.h"
#include "Randomized.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <fstream>

using namespace randomized;
using namespace std;

const size_t seed=1001;

const string filename{"CPPQEDutils_testing_Randomized.d"};

int main(int, char**)
{
  DArray<1> array(10), arrayInterrupted(10), firstHalf(arrayInterrupted(blitz::Range(0,4))), secondHalf(arrayInterrupted(blitz::Range(5,9)));
  
  fillWithRandom(array,seed);
  
  Randomized::Ptr ranFirst(MakerGSL()(seed)), ranSecond(MakerGSL()(seed));
  
  fillWithRandom(firstHalf,ranFirst);
  
  {
    ofstream ofs(filename);
    boost::archive::binary_oarchive ar(ofs);
    ar & *ranFirst;
  }

  // ... some time later restore the other Randomized object instance to its state
  {
    ifstream ifs(filename);
    boost::archive::binary_iarchive ar(ifs);
    ar & *ranSecond;
  }

  fillWithRandom(secondHalf,ranSecond);
  
  return !all(array==arrayInterrupted) ;
  
}
