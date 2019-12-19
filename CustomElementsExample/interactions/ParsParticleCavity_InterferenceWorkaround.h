// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CUSTOMELEMENTSEXAMPLE_INTERACTIONS_PARSPARTICLECAVITY_INTERFERENCEWORKAROUND_H_INCLUDED
#define   CUSTOMELEMENTSEXAMPLE_INTERACTIONS_PARSPARTICLECAVITY_INTERFERENCEWORKAROUND_H_INCLUDED

#include "ModeFunction.h"

#include "ParsFwd.h"

namespace particlecavity_interferenceworkaround {

struct ParsInterference {
  
  ModeFunctionType& modeInterference;
  size_t& kInterference;
  double& uInterference;

  ParsInterference(parameters::ParameterTable&, const std::string& ="");

};

} // particlecavity_interferenceworkaround

#endif // CUSTOMELEMENTSEXAMPLE_INTERACTIONS_PARSPARTICLECAVITY_INTERFERENCEWORKAROUND_H_INCLUDED
