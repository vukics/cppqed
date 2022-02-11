// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "NegPT.h"

#define BOOST_TEST_MODULE NegativityOfThePartialTranspose test
#include <boost/test/unit_test.hpp>

using namespace quantumdata;

typedef details::ExtendV<5,tmptools::Vector<3,1> >::type ExtendedV;


namespace {

using mpl::at_c;
using mpl::equal;

BOOST_STATIC_ASSERT((at_c<ExtendedV,0>::type::value==0));
BOOST_STATIC_ASSERT((at_c<ExtendedV,1>::type::value==6));
BOOST_STATIC_ASSERT((at_c<ExtendedV,2>::type::value==2));
BOOST_STATIC_ASSERT((at_c<ExtendedV,3>::type::value==8));
BOOST_STATIC_ASSERT((at_c<ExtendedV,4>::type::value==4));
BOOST_STATIC_ASSERT((at_c<ExtendedV,5>::type::value==5));
BOOST_STATIC_ASSERT((at_c<ExtendedV,6>::type::value==1));
BOOST_STATIC_ASSERT((at_c<ExtendedV,7>::type::value==7));
BOOST_STATIC_ASSERT((at_c<ExtendedV,8>::type::value==3));
BOOST_STATIC_ASSERT((at_c<ExtendedV,9>::type::value==9));

BOOST_STATIC_ASSERT((equal<ExtendedV,     mpl::vector_c<int,0,6,2,8,4,5,1,7,3,9> >::value));
BOOST_STATIC_ASSERT((equal<ExtendedV,tmptools::Vector  <    0,6,2,8,4,5,1,7,3,9> >::value));

}

BOOST_AUTO_TEST_CASE( Dummy )
{
}

