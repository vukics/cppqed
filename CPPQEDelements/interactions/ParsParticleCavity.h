// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARSPARTICLECAVITY_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARSPARTICLECAVITY_H_INCLUDED

#include "ModeFunction.h"

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

#endif // CPPQEDELEMENTS_INTERACTIONS_PARSPARTICLECAVITY_H_INCLUDED
