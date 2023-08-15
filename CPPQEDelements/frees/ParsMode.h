// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED

// #include "ParsFwd.h"

#include "ComplexExtensions.h"


namespace mode {

struct Pars
{
  size_t &cutoff, &minitFock;
  dcomp& minit;
  size_t &streamLevel;
  double &delta, &omegaKerr, &omegaKerrAlter;

  Pars(parameters::Table&, const std::string& ="");

};


struct ParsDriven : virtual Pars 
{
  dcomp& eta;

  ParsDriven(parameters::Table&, const std::string& ="");

};


struct ParsDissipative : virtual Pars
{
  double &kappa, &nTh;

  ParsDissipative(parameters::Table&, const std::string& ="");

};


struct ParsDrivenDissipative : ParsDriven, ParsDissipative
{
  ParsDrivenDissipative(parameters::Table&, const std::string& ="");
};


} // mode

#endif // CPPQEDELEMENTS_FREES_PARSMODE_H_INCLUDED
