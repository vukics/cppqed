// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_COUPLEDMODESFWD_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_COUPLEDMODESFWD_H_INCLUDED

namespace coupledmodes {

enum Coupling {
  CM_NX,
  CM_XX,
  CM_XX_RWA
};
  
template<bool IS_HA, Coupling C=CM_NX>
class Base;

} // coupledmodes


#endif // CPPQEDELEMENTS_INTERACTIONS_COUPLEDMODESFWD_H_INCLUDED
