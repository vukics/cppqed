// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzTinyExtensions.h"

#include <blitz/tinyvec2io.cc>

#define BOOST_TEST_MODULE BlitzTinyExtensions test
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace blitzplusplus;


template <int I> using IntTiny=blitz::TinyVector<int,I>;


BOOST_AUTO_TEST_CASE( BlitzTinyExtensions )
{
  {
    IntTiny<3> v1{0,1,2};
    IntTiny<4> v2{3,4,5,6};
    
    BOOST_CHECK( blitz::all(concatenateTinies(v1,v2)==IntTiny<7>{0,1,2,3,4,5,6}) ) ;
    BOOST_CHECK( blitz::all(concatenateTinies(v2,v1)==IntTiny<7>{3,4,5,6,0,1,2}) ) ;
  }
  {
    IntTiny<6> v{3,4,2,3,4,2};

    BOOST_CHECK( blitz::all(halfCutTiny(v)==IntTiny<3>{3,4,2}) ) ;
  }
#ifndef   NDEBUG
  IntTiny<6> v{3,4,2,1,4,2};
  
  BOOST_CHECK_THROW( halfCutTiny(v), std::invalid_argument ) ;
#endif // NDEBUG

}
