// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "TMP_Tools.h"

#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>

#include <boost/type_traits/is_base_of.hpp>


using namespace tmptools;


static_assert( Power_v<3,4> ==81 , "Power error" );

static_assert( Power_v<2,18> ==262144 , "Power error" );

// numerical_contains for ...
static_assert( numerical_contains_c<
                                    Range<10,1>,
                                    10ul // ... another type
                                    >::value , "Numerical containment error" );

// 11 should not be in the range 1..10
static_assert( !
               numerical_contains_c<
                                    Range<10,1>,
                                    11
                                    >::value , "Numerical containment error" );

static_assert( numerical_contains_c<
                                    Vector<13,1,87,23,198,2834,21>,
                                    21ul
                                    >::value , "Numerical containment error" );

static_assert( !
               numerical_contains_c<
                                    Vector<13,1,87,23,198,2834,21>,
                                    22u
                                    >::value , "Numerical containment error" );

typedef Vector<87,28,93,1,23,25,97,345,6> V;

static_assert( boost::is_base_of<boost::mpl::vector_c<int,87,28,93,1,23,25,97,345,6>, V >::value , "Is-base-of error" );

static_assert( boost::mpl::size<V>::value==9 , "Size error" );


static_assert( boost::mpl::at_c<V,2>::type::value==93 , "Subscription error" );


typedef pair_c<23,42> P_23_42;


static_assert( P_23_42::first==23 && P_23_42::second==42 , "Pair containment error" );



template struct tmptools::pair_c<23,23,false>;

template struct P_23_42::SanityCheck<21,45>;
