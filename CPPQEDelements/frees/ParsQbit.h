// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED

#include "Qbit_Fwd.h"

#include "ParsFwd.h"

#include "ComplexExtensions.h"


namespace qbit {

struct Pars
{
  dcomp &qbitInit;
  double &delta;

  Pars(parameters::ParameterTable&, const std::string& ="");

};


struct ParsPumped : Pars
{
  dcomp& eta;

  ParsPumped(parameters::ParameterTable&, const std::string& ="");

};


template <typename BASE>
struct ParsLossy : BASE
{
  double &gamma;

  ParsLossy(parameters::ParameterTable&, const std::string& ="");

};


template <typename BASE>
struct ParsLossyPhaseNoise : BASE
{
  double &gamma_parallel;
  
  ParsLossyPhaseNoise(parameters::ParameterTable&, const std::string& ="");
  
};


} // qbit

#endif // CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED
