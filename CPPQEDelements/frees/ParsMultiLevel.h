// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED

#include "MultiLevel_Fwd.h"

#include "ParsFwd.h"

#include "BlitzTiny.h"


namespace multilevel {

/// Type for storing complex level energies (the \f$z_i\f$s \ref multilevelelements "here") \tparam NL number of levels
template<int NL> using ComplexPerLevel = blitz::TinyVector<dcomp,NL>;

/// Type for storing level energies (the \f$\delta_i\f$s \ref multilevelelements "here") \tparam NL number of levels
template<int NL> using RealPerLevel = blitz::TinyVector<double,NL>;


template<int NL, typename VP, typename VL>
struct ParsPumpedLossy
{
  RealPerLevel<NL>& deltas;
  VP& etas;
  VL& gammas;
  double& gamma_parallel;

  ParsPumpedLossy(parameters::ParameterTable&, const std::string& ="");

};


} // multilevel


#endif // CPPQEDELEMENTS_FREES_PARSMULTILEVEL_H_INCLUDED
