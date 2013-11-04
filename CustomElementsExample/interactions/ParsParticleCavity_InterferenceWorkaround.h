// -*- C++ -*-
#ifndef   _PARS_PARTICLE_CAVITY_INTERFERENCE_WORKAROUND_H
#define   _PARS_PARTICLE_CAVITY_INTERFERENCE_WORKAROUND_H

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

#endif // _PARS_PARTICLE_CAVITY_INTERFERENCE_WORKAROUND_H