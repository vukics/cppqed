// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Convenience header for straightforward usage of the trajectory::Simulated class for classical simulations}
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED
#define CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED

#include "BlitzArrayTraits.h"
#include "EvolvedGSL.tcc"
#include "Simulated.tcc"

#include "Pars.tcc"

#include <boost/bind.hpp>

using parameters::ParameterTable;
using parameters::update        ;
using parameters::NamedException;

using trajectory::Simulated;
using trajectory::Pars     ;

#endif // CPPQEDCORE_UTILS_SIMULATED_H_INCLUDED
