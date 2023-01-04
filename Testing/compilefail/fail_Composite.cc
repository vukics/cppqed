// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Composite.h"

template class Composite<composite::result_of::make_list<Act<0,1>,Act<0,2>,Act<0,3>,Act<3,2,5,1> >::type>;
