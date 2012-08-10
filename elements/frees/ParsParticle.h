// -*- C++ -*-
#ifndef   _PARS_PARTICLE_H
#define   _PARS_PARTICLE_H

#include "Particle_Fwd.h"

#include "ModeFunction.h"

#include "ParsFwd.h"


namespace particle {

struct Pars
{

  double &omrec;
  size_t &fin;
  InitialCondition &init;
  int &hoInitn;

  Pars(parameters::ParameterTable&, const std::string& ="");

  virtual ~Pars() {}

};

struct ParsPumped : Pars
{
  double &vClass;
  size_t &kPart;
  ModeFunctionType &modePart;

  ParsPumped(parameters::ParameterTable&, const std::string& ="");

};

} // particle

#endif // _PARS_PARTICLE_H
