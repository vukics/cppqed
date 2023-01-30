// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Archive.h"
#include "Random.h"

#define BOOST_TEST_MODULE RandomEngineSerialization test
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <span>
#include <vector>

using namespace std;

const size_t seed=1001;

template<typename RandomEngine>
auto testcase()
{
  vector<double>
    array(10),
    arrayInterrupted(10);
    
  span
    firstHalf{arrayInterrupted.begin(),arrayInterrupted.begin()+5},
    secondHalf{arrayInterrupted.begin()+5,arrayInterrupted.end()};
  
  randomutils::fill<uniform_real_distribution<double>,RandomEngine>(array,seed);
  
  cerr<<cppqedutils::toStringJSON(array)<<endl;
  
  RandomEngine reFirst(seed), reSecond(seed);
  
  randomutils::fill<uniform_real_distribution<double>>(firstHalf,reFirst);
  
  const std::string filename{"CPPQEDutils_testing_"+randomutils::EngineID_v<RandomEngine>+".d"};
  
  BOOST_TEST_CHECKPOINT( "PRE_OUTPUT" );
  
  {
    ofstream ofs(filename);
    cppqedutils::oarchive ar(ofs, ios_base::trunc | ios_base::binary);
    ar & reFirst;
  }

  BOOST_TEST_CHECKPOINT( "POST_OUTPUT" );

  // ... some time later restore the other RandomEngine instance to its state
  {
    ifstream ifs(filename);
    cppqedutils::iarchive ar(ifs, ios_base::binary);
    ar & reSecond;
  }

  BOOST_TEST_CHECKPOINT( "POST_INPUT" );

  randomutils::fill<uniform_real_distribution<double>>(secondHalf,reSecond);

  cerr<<cppqedutils::toStringJSON(arrayInterrupted)<<endl;

  return array==arrayInterrupted ;
}


// BOOST_AUTO_TEST_CASE( GSL_RandomEngineTest ) { BOOST_CHECK(testcase<randomutils::GSL_Engine>()); }

BOOST_AUTO_TEST_CASE( STD_MersenneTwisterTest ) { BOOST_CHECK(testcase<mt19937_64>()); }

BOOST_AUTO_TEST_CASE( PCG64_Test ) { BOOST_CHECK(testcase<pcg64>()); }

BOOST_AUTO_TEST_CASE( Xoshiro256pp_Test ) { BOOST_CHECK(testcase<XoshiroCpp::Xoshiro256PlusPlus>()); }
