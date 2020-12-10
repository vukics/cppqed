// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_H_INCLUDED

#include "Pars.h"

#include <boost/fusion/sequence/io.hpp>
#include <boost/fusion/include/io.hpp>


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
