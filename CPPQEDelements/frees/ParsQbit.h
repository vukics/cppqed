// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
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


struct ParsPumped : virtual Pars 
{
  dcomp& eta;

  ParsPumped(parameters::ParameterTable&, const std::string& ="");

};


struct ParsLossy : virtual Pars
{
  double &gamma;

  ParsLossy(parameters::ParameterTable&, const std::string& ="");

};


struct ParsPumpedLossy : ParsPumped, ParsLossy
{
  ParsPumpedLossy(parameters::ParameterTable&, const std::string& ="");
};


} // qbit

#endif // CPPQEDELEMENTS_FREES_PARSQBIT_H_INCLUDED
