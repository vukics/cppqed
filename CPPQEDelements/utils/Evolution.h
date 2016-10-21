// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Convenience header collecting to one place most of the components of C++QED core needed for writing scripts (basically, the core part of the Level-1 interface)}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED

#include "BlitzArrayTraits.h"

#include "Evolution.tcc"
#include "QM_Picture.h"
#include "Tridiagonal.tcc"

#include "EvolvedGSL.tcc"
#include "Pars.tcc"

#include "component_versions.h"


/// Introduces ParameterTable into the global namespace to break ambiguity between update and parameters::update
class ParameterTable : public parameters::ParameterTable {};


/// Convenience version of parameters::update that includes the highest-level version information
/** \note This cannot be put into an implementation file because then an incorrect `component_versions.h` file would be picked up */
void update(ParameterTable& p, int argc, char* argv[], const std::string& prefix="--")
{
  updateVersionstring(cppqed_component_versions());
  parameters::update(p,argc,argv,prefix);
}

/// Convenience version of picture::updateWithPicture that includes the highest-level version information
QM_Picture& updateWithPicture(ParameterTable& p, int argc, char* argv[], const std::string& prefix="--")
{
  updateVersionstring(cppqed_component_versions());
  return picture::updateWithPicture(p,argc,argv,prefix);
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED

