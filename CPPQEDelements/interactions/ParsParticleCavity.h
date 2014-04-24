// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARSPARTICLECAVITY_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARSPARTICLECAVITY_H_INCLUDED

#include "ModeFunction.h"

#include "ParsFwd.h"


namespace particlecavity {

struct ParsGenericPump {

  double &etaeff;

protected:
  ParsGenericPump(parameters::ParameterTable&, const std::string& ="");
};

struct ParsOrthogonal {

  double &uNot;

  ParsOrthogonal(parameters::ParameterTable&, const std::string& ="");

};

struct ParsAlong : ParsOrthogonal {

  size_t &kCav;
  ModeFunctionType &modeCav;

  ParsAlong     (parameters::ParameterTable&, const std::string& ="");

};

struct ParsOrthogonalGenericPump : ParsOrthogonal, ParsGenericPump {
  ParsOrthogonalGenericPump(parameters::ParameterTable &p, const std::string &mod="")
    : ParsOrthogonal(p,mod), ParsGenericPump(p,mod) {}
};

struct ParsAlongGenericPump : ParsAlong, ParsGenericPump {
  ParsAlongGenericPump(parameters::ParameterTable &p, const std::string &mod="")
    : ParsAlong(p,mod), ParsGenericPump(p,mod) {}
};

} // particlecavity

#endif // CPPQEDELEMENTS_INTERACTIONS_PARSPARTICLECAVITY_H_INCLUDED
