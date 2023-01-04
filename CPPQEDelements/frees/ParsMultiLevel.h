// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED

#include "Pars.h"

#include "BlitzTiny.h"


namespace multilevel {

/// Type for storing complex level energies (the \f$z_i\f$s \ref multilevelelements "here") \tparam NL number of levels
template<int NL> using ComplexPerLevel = blitz::TinyVector<dcomp,NL>;

/// Type for storing level energies (the \f$\delta_i\f$s \ref multilevelelements "here") \tparam NL number of levels
template<int NL> using RealPerLevel = blitz::TinyVector<double,NL>;


template<int NL, typename VP, typename VL>
struct ParsDrivenDissipative
{
  RealPerLevel<NL>& deltas;
  VP& etas;
  VL& gammas;
  double& gamma_parallel;

  RealPerLevel<NL>& gamma_parallel_vector;
  
  ParsDrivenDissipative(parameters::Table& p, const std::string& mod="")
    : deltas(p.addTitle("DrivenDissipativeMultiLevel",mod).add("deltas",mod,"MultiLevel detunings vector",RealPerLevel<NL>{0.})),
      etas(p.add("etas",mod,"MultiLevel pumps vector",VP())),
      gammas(p.add("gammas",mod,"MultiLevel decays vector",VL())),
      gamma_parallel(p.add("gamma_parallel",mod,"Phase flip rate",0.)),
      gamma_parallel_vector(p.add("gamma_parallel_vector",mod,"Phase flip rate when it’s different for the individual levels",RealPerLevel<NL>{0.}))
  {
    deltas=0;
  }

};


} // multilevel


#endif // CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED
