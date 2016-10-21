// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED

#include "Mode_Fwd.h"

#include "ParsFwd.h"

#include "ComplexExtensions.h"


namespace mode {

struct Pars
{
  size_t &cutoff, &minitFock;
  dcomp& minit;
  size_t &displayLevel;
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
  double &kappa, &nTh;

  ParsLossy(parameters::ParameterTable&, const std::string& ="");

};


struct ParsPumpedLossy : ParsPumped, ParsLossy
{
  ParsPumpedLossy(parameters::ParameterTable&, const std::string& ="");
};


} // mode

#endif // CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED
