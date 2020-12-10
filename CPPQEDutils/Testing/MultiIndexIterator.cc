// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MultiIndexIterator.h"
#define BOOST_TEST_MODULE MultiIndexIterator test
#include <boost/test/unit_test.hpp>

using namespace cpputils

typedef IdxTiny<3> Idx;

// BOOST_AUTO_TEST_SUITE( test_suite )

BOOST_AUTO_TEST_CASE( iterator )
{
  Idx lbound(-2,0,1), ubound(-1,1,3);

  MultiIndexIterator<3> 
    begin(lbound,ubound,std::false_type()),
    end  (lbound,ubound,std:: true_type());

  BOOST_CHECK(all(*begin++==lbound));

  BOOST_CHECK(all(*begin++==Idx(-2,0,2)));
  BOOST_CHECK(all(*begin++==Idx(-2,0,3)));
  BOOST_CHECK(all(*begin++==Idx(-2,1,1)));
  BOOST_CHECK(all(*begin++==Idx(-2,1,2)));
  BOOST_CHECK(all(*begin++==Idx(-2,1,3)));
  BOOST_CHECK(all(*begin++==Idx(-1,0,1)));
  BOOST_CHECK(all(*begin++==Idx(-1,0,2)));
  BOOST_CHECK(all(*begin++==Idx(-1,0,3)));
  BOOST_CHECK(all(*begin++==Idx(-1,1,1)));
  BOOST_CHECK(all(*begin++==Idx(-1,1,2)));
  BOOST_CHECK(all(*begin++==Idx(-1,1,3)));

  BOOST_CHECK(begin++==end);
  
  BOOST_CHECK(all(*end    ==Idx(0,0,1)));  

  BOOST_CHECK(all(*begin  ==Idx(0,0,2)));

}

// BOOST_AUTO_TEST_SUITE_END()


