// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Combinatorics.h"
#include "ContainerIO.h"

#define BOOST_TEST_MODULE CWR test
#include <boost/test/unit_test.hpp>

using namespace cppqedutils;
using namespace std;


BOOST_AUTO_TEST_CASE( CWR )
{
  CWR_Dir<4> dir{5};

  cout<<dir()<<endl;

  const auto c12=dir[12];
  BOOST_CHECK(dir[c12]==12);

  BOOST_CHECK( ( dir[{2,2,0,1}]==43 ) );

  BOOST_CHECK_THROW( ( dir[{1,2,0,1}] ) , CWR_Dir<4>::SubscriptingException );

  BOOST_CHECK_THROW( ( dir[{2,2,1,1}] ) , CWR_Dir<4>::SubscriptingException );

  BOOST_CHECK_THROW( ( dir[{2,3,1}] ) , CWR_Dir<4>::SubscriptingException );

  // const boost::array<size_t,5> idxBoostErr={{2,2,0,1,2}};

}
