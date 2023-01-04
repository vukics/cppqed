// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED

#include "Pars.h"


namespace mljc {

template<typename VC>
struct Pars
{
  VC& gs;

  Pars(parameters::Table& p, const std::string& mod="") : gs(p.addTitle("MultiLevelJaynesCummings",mod).add("gs",mod,"Multi-Level Jaynes-Cummings couplings",VC()))
  {}

};

} // mljc



#endif // CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED
