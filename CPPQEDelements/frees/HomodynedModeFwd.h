// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_HOMODYNEDMODEFWD_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_HOMODYNEDMODEFWD_H_INCLUDED

#include "Mode_Fwd.h"

template<typename A=mode::Averaged>
class HomodynedMode;

namespace mode {

template<typename BASE>
struct ParsHomodyned; // BASE either ParsLossy or ParsPumpedLossy

} // mode

#endif // CPPQEDELEMENTS_FREES_HOMODYNEDMODEFWD_H_INCLUDED
