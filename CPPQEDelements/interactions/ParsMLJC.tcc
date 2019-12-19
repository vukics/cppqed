// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_TCC_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_TCC_INCLUDED

#include "Pars.tcc"

#include <boost/fusion/sequence/io.hpp>
#include <boost/fusion/include/io.hpp>


namespace mljc {


template<typename VC>
Pars<VC>::Pars(parameters::ParameterTable& p, const std::string& mod)
  : gs(p.addTitle("MultiLevelJaynesCummings",mod).addMod("gs",mod,"Multi-Level Jaynes-Cummings couplings",VC()))
{
}


} // mljc


#endif // CPPQEDELEMENTS_INTERACTIONS_PARSMLJC_TCC_INCLUDED
