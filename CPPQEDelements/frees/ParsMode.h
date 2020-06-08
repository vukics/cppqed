// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED

#include "ParsFwd.h"

#include "ComplexExtensions.h"


namespace mode {

struct Pars
{
  size_t &cutoff, &minitFock;
  dcomp& minit;
  size_t &displayLevel;
  double &delta;

  Pars(parameters::Table&, const std::string& ="");

};


struct ParsPumped : virtual Pars 
{
  dcomp& eta;

  ParsPumped(parameters::Table&, const std::string& ="");

};


struct ParsLossy : virtual Pars
{
  double &kappa, &nTh;

  ParsLossy(parameters::Table&, const std::string& ="");

};


struct ParsPumpedLossy : ParsPumped, ParsLossy
{
  ParsPumpedLossy(parameters::Table&, const std::string& ="");
};


} // mode

#endif // CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED
