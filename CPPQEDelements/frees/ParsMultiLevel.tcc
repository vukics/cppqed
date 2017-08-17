// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_TCC_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_TCC_INCLUDED

#include "Pars.tcc"


namespace multilevel {


template<int NL, typename VP, typename VL>
ParsPumpedLossy<NL,VP,VL>::ParsPumpedLossy(parameters::ParameterTable& p, const std::string& mod)
  : deltas(p.addTitle("PumpedLossyMultiLevel",mod).addMod("deltas",mod,"MultiLevel detunings vector",Levels())),
    etas(p.addMod("etas",mod,"MultiLevel pumps vector",VP())),
    gammas(p.addMod("gammas",mod,"MultiLevel decays vector",VL())),
    gamma_parallel(p.addMod("gamma_parallel",mod,"Phase flip rate",0.))
{
  deltas=0;
}


} // multilevel


#endif // CPPQEDELEMENTS_FREES_PARSMULTILEVEL_TCC_INCLUDED
