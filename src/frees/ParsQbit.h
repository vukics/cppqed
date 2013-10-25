// -*- C++ -*-
#ifndef   ELEMENTS_FREES_PARSQBIT_H_INCLUDED
#define   ELEMENTS_FREES_PARSQBIT_H_INCLUDED

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

#endif // ELEMENTS_FREES_PARSQBIT_H_INCLUDED
