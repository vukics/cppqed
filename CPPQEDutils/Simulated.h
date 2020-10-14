// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Convenience header for straightforward usage of the trajectory::Simulated class for classical simulations}
#ifndef CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED
#define CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED

#include "BlitzArray.h"
#include "EvolvedGSL.tcc"
#include "Simulated_.h"

#include "Pars.h"

using ParameterTable=parameters::Table;
using parameters::update        ;
using parameters::NamedException;

using trajectory::Simulated;
using trajectory::Pars     ;

#endif // CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED
