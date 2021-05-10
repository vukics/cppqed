// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "TMP_Tools.h"

#include <boost/type_traits/is_base_of.hpp>


using namespace tmptools;


static_assert( Power_v<3,4> ==81 , "Power error" );

static_assert( Power_v<2,18> ==262144 , "Power error" );

typedef Vector<87,28,93,1,23,25,97,345,6> V;


typedef pair_c<23,42> P_23_42;


static_assert( P_23_42::first==23 && P_23_42::second==42 , "Pair containment error" );



template struct tmptools::pair_c<23,23,false>;

template struct P_23_42::SanityCheck<21,45>;
