// -*- C++ -*-
#ifndef   ELEMENTS_FREES_PARSMODE_H_INCLUDED
#define   ELEMENTS_FREES_PARSMODE_H_INCLUDED

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

#endif // ELEMENTS_FREES_PARSMODE_H_INCLUDED
