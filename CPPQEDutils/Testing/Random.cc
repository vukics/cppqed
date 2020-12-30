// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Archive.h"
#include "BlitzArray.h"
#include "Random.h"

#include <fstream>

using namespace cpputils;
using namespace std;

const size_t seed=1001;

const string filename{"CPPQEDutils_testing_Random.d"};

int main(int, char**)
{
  DArray<1>
    array(10),
    arrayInterrupted(10),
    firstHalf(arrayInterrupted(blitz::Range(0,4))),
    secondHalf(arrayInterrupted(blitz::Range(5,9)));
  
  fillWithRandom<std::uniform_real_distribution<double>,cpputils::GSL_RandomEngine>(array,seed);
  
  std::cerr<<array<<std::endl;
  
  cpputils::GSL_RandomEngine reFirst(seed), reSecond(seed);
  
  fillWithRandom<std::uniform_real_distribution<double>>(firstHalf,reFirst);
  
  std::cerr<<array<<std::endl;
  
  {
    ofstream ofs(filename);
    oarchive ar(ofs);
    ar & reFirst;
  }

  // ... some time later restore the other RandomEngine instance to its state
  {
    ifstream ifs(filename);
    iarchive ar(ifs);
    ar & reSecond;
  }

  fillWithRandom<std::uniform_real_distribution<double>>(secondHalf,reSecond);

  std::cerr<<arrayInterrupted<<std::endl;

  return !all(array==arrayInterrupted) ;
  
}
