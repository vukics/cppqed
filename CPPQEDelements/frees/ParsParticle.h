// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_PARSPARTICLE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARSPARTICLE_H_INCLUDED

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

#endif // CPPQEDELEMENTS_FREES_PARSPARTICLE_H_INCLUDED
