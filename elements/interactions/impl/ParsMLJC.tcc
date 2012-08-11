// -*- C++ -*-
#ifndef   PARS_MULTI_LEVEL_JAYNES_CUMMINGS_IMPL_INCLUDED
#define   PARS_MULTI_LEVEL_JAYNES_CUMMINGS_IMPL_INCLUDED

#include "impl/Pars.tcc"

#include <boost/fusion/sequence/io.hpp>
#include <boost/fusion/include/io.hpp>


namespace mljc {


template<typename VC>
Pars<VC>::Pars(parameters::ParameterTable& p, const std::string& mod)
  : gs(p.addTitle("MultiLevelJaynesCummings",mod).addMod("gs",mod,"Multi-Level Jaynes-Cummings couplings",VC()))
{
}


} // mljc


#endif // PARS_MULTI_LEVEL_JAYNES_CUMMINGS_IMPL_INCLUDED
