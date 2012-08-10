// -*- C++ -*-
#ifndef   _PARS_MODE_H
#define   _PARS_MODE_H

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

#endif // _PARS_MODE_H
