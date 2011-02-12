// -*- C++ -*-
#ifndef   _PARS_PARTICLE_CAVITY_H
#define   _PARS_PARTICLE_CAVITY_H

#include "ModeFunctionTypeFwd.h"

#include "ParsFwd.h"

namespace particlecavity {

struct ParsOrthogonal {
  
  double &uNot;
  
  ParsOrthogonal(parameters::ParameterTable&, const std::string& ="");

};


struct ParsAlong : ParsOrthogonal {
  
  size_t &kCav;
  ModeFunctionType &modeCav;
  
  ParsAlong     (parameters::ParameterTable&, const std::string& ="");

};


} // particlecavity

#endif // _PARS_PARTICLE_CAVITY_H
