// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Combinatorics.h"

#define BOOST_TEST_MODULE CWR test
#include <boost/test/unit_test.hpp>

#include "BlitzTiny.h"

#include <boost/array.hpp>

using namespace cpputils;
using namespace std;


BOOST_AUTO_TEST_CASE( CWR )
{
  CWR_Dir dir(4,5);

  cout<<dir()<<endl;

  const CWR_Dir::Configuration c12=dir[12];
  BOOST_CHECK(dir[c12]==12);

  const IdxTiny<4> idxTiny(2,2,0,1);
  
  const boost::array<size_t,4> idxBoost{{2,2,0,1}};
  BOOST_CHECK(dir[idxBoost]==12);

  // const boost::array<size_t,5> idxBoostErr={{2,2,0,1,2}};

  BOOST_CHECK(equal(c12.begin(),c12.end(),idxTiny.begin()));

  std::vector<size_t> idxSTL(c12.begin(),c12.end());
  BOOST_CHECK(dir[idxSTL]==12);

  const CWR_Dir::Configuration idxBlitz(c12.copy());
  BOOST_CHECK(dir[idxBlitz]==12);

}
