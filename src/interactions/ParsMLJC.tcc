// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_IMPL_PARSMLJC_TCC_INCLUDED
#define   ELEMENTS_INTERACTIONS_IMPL_PARSMLJC_TCC_INCLUDED

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


#endif // ELEMENTS_INTERACTIONS_IMPL_PARSMLJC_TCC_INCLUDED